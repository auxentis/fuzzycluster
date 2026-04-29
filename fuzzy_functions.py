import re
import numpy as np
from operator import itemgetter
import os,sys
import math
import matplotlib.pyplot as plt
import pandas as pd
import csv
from termcolor import colored
from Bio.SeqIO.FastaIO import FastaIterator
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import io
from wtforms import ValidationError

def pattern_element_length(min_len=3, max_len=30, patternsnum=3):
    """Validator that ensures each element separated by '|' is <= max_len and > min_len."""
    def _validator(form, field):
        raw = field.data.split('|')
        if any(p.strip() == "" for p in raw):
             raise ValidationError("Empty patterns are not allowed. Avoid leading or trailing '|'.")

        patterns = [p.strip() for p in raw]

        if len(patterns) < 1 or len(patterns) > patternsnum:
                raise ValidationError(
                f"Up to {patternsnum} different patterns can be analyzed at present, entered are {len(patterns)} patterns"
                )   
        else: 
            for p in patterns:
                if len(p) > max_len or len(p) < min_len:
                    raise ValidationError(
                        f"Each pattern must be minimum {min_len} and at most {max_len} characters "
                        f"'{p}' is {len(p)} characters long."
                )     
            
    return _validator


def fasta_reader(infile):
    from Bio.SeqIO.FastaIO import FastaIterator
    with open(infile) as handle:
        for record in FastaIterator(handle):
            yield record



def search_functionGC(myseq, pattern):
    ind_spans = [match.span() for match in re.finditer(pattern, myseq)]

    ind_start=[span[0] for span in ind_spans]
    ind_end=[span[1] for span in ind_spans]

    spans_length = np.subtract(ind_end, ind_start)
    if len(spans_length) > 0:
        spans_median = np.median(spans_length)
    else:
        spans_median = 0
    
    return ind_spans, ind_start, ind_end, spans_median

def highlight_functionGC(myseq, pattern_span_open, span_close, ind_spans):

    if len(ind_spans)>0: 

        highlighted_list=list(myseq)
        
        for start, end in ind_spans:
            start = start
            end = end-1
            highlighted_list[start]=pattern_span_open+highlighted_list[start]
            highlighted_list[end]=highlighted_list[end] +span_close
        highlighted=''.join(highlighted_list)


    else:
        ind_spans ='None'
        highlighted = 'None'

    return highlighted


def loop_function(myseq, start, end):
    """
    Determines the loop lengths between repeats

    """

    ind_start = np.array(start)
    ind_end = np.array(end)
    sequence_length = len(myseq)

    if len(ind_start)>0:
        loop_list = np.subtract(ind_start[1:], ind_end[:-1])
        start_loop = ind_start[0]
        end_loop = sequence_length -ind_end[-1]

        list_zeros = np.zeros(len(loop_list), dtype=int)
        loop_list_zeroes = (np.dstack((loop_list, list_zeros)))[0].flatten()
        loop_list_zeroes_stack = np.insert((np.dstack((loop_list_zeroes, loop_list_zeroes)))[0].flatten(), [0,0,0,0], [start_loop, start_loop, 0,0])
        loop_list_end_loops = np.append(loop_list_zeroes_stack, (end_loop, end_loop))
        
        loop_list_left_tail = np.insert(loop_list, 0, start_loop)
        loop_list_arr = np.append(loop_list_left_tail, end_loop)
        loop_median = np.median(loop_list_arr)
        loop_ind = (np.dstack((ind_start, ind_end)))[0]
        loop_ind_din = np.insert((np.dstack((loop_ind, loop_ind))).flatten(), 0,0)
        loop_ind_din_appended = np.append(loop_ind_din, sequence_length)
        
    else:
        loop_ind_din_appended = None
        loop_list_end_loops = None
        loop_list_arr = None
        loop_median = sequence_length

    return loop_ind_din_appended, loop_list_end_loops, loop_list_arr, loop_median


def merge_intervals(intervals):
    """
    Merging of Tel intervals with loops to get Fussy Clusters

    """
    starts = intervals[:,0]
    ends = np.maximum.accumulate(intervals[:,1])
    valid = np.zeros(len(intervals) + 1, dtype=bool)
    valid[0] = True
    valid[-1] = True
    valid[1:-1] = starts[1:] >= ends[:-1]
    fu_cluster = np.vstack((starts[:][valid[:-1]], ends[:][valid[1:]])).T 
    valid11 = valid[1:-1].tolist()
    return [fu_cluster, valid11]

def cluster_functionGC(myseq, name, clusterscore, start, end, loop_list_arr, input_looplength, loop_median, spans_median, input_ssr, input_repeats, tandem_select):
    """
    Calculates clusters depending on provided parameters
    Calculates Cluster Score, SSR, TCS
    Makes filtration
    Outputs clusters dataframe and loops

    """
    if clusterscore == 'undefined':
        clusterscore = 0
    else:
        clusterscore = clusterscore
    
    if tandem_select == True:
        input_repeats_number = 2
    else:
        if input_repeats == 'any':
            input_repeats_number = 2    
        else:
            input_repeats_number = int(input_repeats)

    if input_ssr == 'undefined':
        input_ssr_ratio = float(1.5)
    else:
        input_ssr_ratio = input_ssr

    if input_looplength == 'estimate':
        loop_estimate = round(loop_median)
    else:
        loop_estimate = input_looplength

    ind_start = np.array(start)
    ind_end = np.array(end)
    sequence_length = len(myseq)
    loop_list_complete = loop_list_arr
    
    if len(ind_start)>0:
        interloops = loop_list_complete[1:-1]
        loop_condition = interloops<= loop_estimate
        ind_end_condition = np.where(loop_condition, ind_end[1:], ind_end[:-1])
        ind_end_adj = np.append(ind_end_condition, ind_end[-1])
        zipped_list = (np.dstack((ind_start, ind_end_adj)))[0]
        info_loop_split = np.split(interloops, np.where(interloops > loop_estimate)[0])
        total_loops_len = [np.sum (i[i<=loop_estimate]) for i in info_loop_split]

        #new res_leni
        res_leni = [i.size for i in info_loop_split]
        res_leni[0] = res_leni[0]+1
    

        if len(interloops) >= 1:
            fused_cluster_join = merge_intervals(zipped_list)
            fused_cluster = fused_cluster_join[0]
            valid_list = fused_cluster_join[1]
    
            #Determining clusters length
            a = fused_cluster[:,1]
            b = fused_cluster[:,0]
            fused_interval = np.subtract(a, b)
            
            if len(res_leni) >=1:
        
            # #Counting Score - a parameter which reflects the number of Tel intervals to the length of Tel Fused Cluster 
            # (Score = (n(Tel)^2/l(Fused Cluster))*sqrt(l))
                fused_score = np.divide(np.power(res_leni, 2), np.sqrt(fused_interval))
                theoretical_cluster_length = np.add(np.subtract(fused_interval, total_loops_len), np.multiply(np.subtract(res_leni, 1), loop_estimate))
                theoretical_score = np.divide(np.power(res_leni, 2), np.sqrt(theoretical_cluster_length))
                ssr = np.divide(theoretical_score, fused_score)
                fused_score = np.around(fused_score, 2)
                ssr = np.around(ssr, 3)
            else:
                pass
            
            cluster_df_sorted=pd.DataFrame({   
               "Sequence name": [name]*len(fused_cluster),
               "Cluster Start": fused_cluster[:,0],
               "Cluster End" : fused_cluster[:,1],
               "Cluster length" : fused_interval,
               "Number of repeats": res_leni,
               "Total Loops length": total_loops_len,
               "Score": fused_score,
               "Input loop": loop_estimate,
               "SSR": ssr})

            
            # Filtration
            score_limit = float(clusterscore)
            cluster_df_select = cluster_df_sorted[(cluster_df_sorted["Score"] >= score_limit) & 
                                                  (cluster_df_sorted["SSR"] <= input_ssr_ratio) & 
                                                  (cluster_df_sorted["Number of repeats"] >= input_repeats_number)]
            
            if cluster_df_select.empty:
                cluster_df_select = pd.DataFrame(np.zeros((1, 9), dtype=int), columns=['Sequence name', 'Cluster Start', 'Cluster End', 'Cluster length', 'Number of repeats', 'Total Loops length', 'Score', 'Input loop', 'SSR'])
                cluster_df_select['Sequence name']=name
            else:
                cluster_df_select = cluster_df_select
        
        else:

            cluster_df_select = pd.DataFrame(np.zeros((1, 9), dtype=int), columns=['Sequence name', 'Cluster Start', 'Cluster End', 'Cluster length', 'Number of repeats', 'Total Loops length', 'Score', 'Input loop', 'SSR'])
            cluster_df_select['Sequence name']=name
    
    else:
        cluster_df_select = pd.DataFrame(np.zeros((1, 9), dtype=int), columns=['Sequence name', 'Cluster Start', 'Cluster End', 'Cluster length', 'Number of repeats', 'Total Loops length', 'Score', 'Input loop', 'SSR'])
        cluster_df_select['Sequence name']=name
    
    blankIndex=[''] * len(cluster_df_select)
    cluster_df_select.index=blankIndex

    return cluster_df_select, loop_estimate


def cluster_loopsGC(myseq, cluster_array):

    sequence_length = len(myseq)
    cluster_array = cluster_array.to_numpy()
    cluster_start = cluster_array[:,0]
    cluster_end = cluster_array[:,1]
    cluster_length = cluster_array[:,2]

    zeros_array = np.zeros(len(cluster_length), dtype=int)
    cluster_zeropadded = (np.dstack((cluster_length, zeros_array)))[0]
    cluster_flat = cluster_zeropadded.flatten()
    cluster_indexes_stack = (np.dstack((cluster_start, cluster_end)))[0]
    cluster_indexes = cluster_indexes_stack.flatten()
    cluster_indexes_list = cluster_indexes_stack.tolist()
    cluster_list = np.insert((np.dstack((cluster_flat, cluster_flat)))[0].flatten(),[0, 0], [0, 0])
    
    cluster_indzip = np.insert((np.dstack((cluster_indexes, cluster_indexes)))[0].flatten(), 0, 0)
    cluster_ind = np.append(cluster_indzip, sequence_length)
    
    return cluster_list, cluster_ind, cluster_indexes_list


def highlighted_clusterGC(myseq, cluster_indexes, cluster_list, cluster_span, cluster_span_open, span_close, pattern):

    if sum(cluster_list) >0:

        highlighted_list=list(myseq)
    
        for start, end in cluster_indexes:
            start = start
            end = end-1
            highlighted_list[start]=cluster_span_open+highlighted_list[start]
            highlighted_list[end]=highlighted_list[end] +span_close

        highlighted_cluster=''.join(highlighted_list)
        highlighted_cluster_repeat = re.sub(pattern, cluster_span, highlighted_cluster)

    
    else:
        highlighted_cluster_repeat = re.sub(pattern, cluster_span, myseq)

    return highlighted_cluster_repeat

    

def cluster_plot(myseq, clusterG_list, clusterG_ind, clusterC_list, clusterC_ind, description,save_dir, selected_option):
    """
    Makes cluster plots
    """
    if selected_option == 'TTAGGG':
        labelF = 'ClusterG'
        labelR = 'ClusterC'
    
    elif selected_option == 'FuzzyTel':
        labelF = 'ClusterG'
        labelR = 'ClusterC'

    elif selected_option == 'Custom pattern':
        labelF = 'ClusterF'
        labelR = 'ClusterR'
    
    xlen = len(myseq)
    
    set_fontsize = 8
        
    def fig1():
        if sum(clusterG_list) > 0:
            if sum(clusterC_list) > 0:
                fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(clusterG_ind, clusterG_list, label=labelF, color=(0,0,1))
                ax1.legend()
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Cluster length")
                            
                ax2.plot(clusterC_ind, clusterC_list, label=labelR, color=(0,1,0))
                ax2.legend()
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Cluster length")
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig1.suptitle('Map of Repeat clusters')
                fig1.tight_layout()
                return fig1
            else:
                fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(clusterG_ind, clusterG_list, label=labelF, color=(0,0,1))
                ax1.legend()
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Cluster length")
                            
                ax2.plot(0, 0, label=labelR, color=(0,1,0))
                ax2.legend()
                ax2.set_ylim(0, 100)
                ax2.set_xlim(0, xlen)
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Loop length")
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig1.suptitle('Map of Repeat clusters')
                fig1.tight_layout()
                return fig1

        elif sum(clusterC_list) > 0:    
            fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Cluster length")
                            
            ax2.plot(clusterC_ind, clusterC_list, label=labelR, color=(0,1,0))
            ax2.legend()
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Cluster length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig1.suptitle('Map of Repeat clusters')
            fig1.tight_layout()
            return fig1
        else:
            fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Cluster length")
                            
            ax2.plot(0, 0, label=labelR, color=(0,1,0))
            ax2.legend()
            ax2.set_ylim(0, 100)
            ax2.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Cluster length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig1.suptitle('Map of Repeat clusters')
            fig1.tight_layout()
            return fig1

    figure1 = fig1()
    figure1.savefig(os.path.join(save_dir, 'cluster_plot.png'), dpi=100)
    plt.close(figure1)

    # Derive session folder name
    session_folder = os.path.basename(save_dir)
    # Create URLs compounds that Flask can serve
    plot3_url = f"{session_folder}/cluster_plot.png"
    return plot3_url  


def bin_range(loop_list_arrG, loop_list_arrC):
    """
    Calculates the optimal bin range for a given list of coordinates
    for complementary strands, presented as subplots. Thus, that bins are equal for the subplots.

    Args:
    -----
    loop_list_arrG, loop_list_arrC : `array <int>`
        Numpy 1-d array with genomic coordinates
    Returns:
    --------
    bins : `array <int>`
        The list of bins, ready for plotting
    """
    if loop_list_arrG is None:
        maxG = 0
        minG = 0
        lenG = 0
        sdevG = 0
        percent99G = 0
    else:
        maxG = np.max(loop_list_arrG)
        minG = np.min(loop_list_arrG)
        lenG = len(loop_list_arrG)
        sdevG = np.std(loop_list_arrG)
        percent99G = np.percentile(loop_list_arrG, 99)   

    if loop_list_arrC is None:
        maxC = 0
        minC = 0
        lenC = 0
        sdevC = 0
        percent99C = 0
    else:  
        maxC = np.max(loop_list_arrC)
        minC = np.min(loop_list_arrC)
        lenC = len(loop_list_arrC)
        sdevC = np.std(loop_list_arrC)
        percent99C = np.percentile(loop_list_arrC, 99)
    maxGC = max(maxG, maxC)
    maxdev = max(sdevG, sdevC)
    maxpercent = max(percent99G, percent99C)
    lenGC = max(lenG, lenC)   
    if maxdev+maxpercent >= maxGC:
        maxpercent_adj = maxdev+maxpercent
    else:
        maxpercent_adj = maxpercent
    
    if maxdev+maxpercent == 0:
        maxpercent_adj = lenGC

    if lenGC <= 4:
        max_bin_count = 10
    else:
        max_bin_count = min(100, lenGC)

    step=math.ceil(maxpercent_adj/max_bin_count)
    maxpercent_plus = maxpercent_adj+step
    bins_range = np.arange(0,maxpercent_plus,step)
    bins_range_plus = np.append(bins_range, bins_range[-1]+step)
    
    return bins_range_plus


def save_image(file, fig_list):

    pdf = PdfPages(file)

    for fig in fig_list: 
        fig.savefig(pdf, format='pdf') 
          
    pdf.close()  


def image_functionGC(myseq, loop_indG, loop_listG, loop_list_arrG, loop_indC, loop_listC, loop_list_arrC, description, loop_medianG, loop_medianC,save_dir, selected_option):
    """
    Draws loops map and hist map

    """
    if selected_option == 'TTAGGG':
        labelF = 'LoopG'
        labelR = 'LoopC'
    
    elif selected_option == 'FuzzyTel':
        labelF = 'LoopG'
        labelR = 'LoopC'

    elif selected_option == 'Custom pattern':
        labelF = 'LoopF'
        labelR = 'LoopR'
    
    set_fontsize = 8

    xlen = len(myseq)

    # degenerate_patterns='[NRYBDKMHVSW]+'
    degenerate_patterns = '[^ATGC]+'

    unknown_pattern=re.compile(degenerate_patterns,re.IGNORECASE)
    unknown_ind=[( x.start(), x.end() ) for x in re.finditer(unknown_pattern,myseq)]


    def fig1():

        if loop_indG is not None:
            if loop_indC is not None:
                maxG_y = max(loop_listG)
                maxC_y = max(loop_listC)

                fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(loop_indG, loop_listG, label=labelF, color=(0,0,1))
                unknown_level=float(maxG_y)*1.1
                first_unknown=True
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1  
                    ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax1.legend()
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Loop length")
                        
                ax2.plot(loop_indC, loop_listC, label=labelR, color=(0,1,0))
                unknown_level=float(maxC_y)*1.1
                first_unknown=True
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1  
                    ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax1.autoscale(enable=True, axis='x', tight=False)
                ax2.autoscale(enable=True, axis='x', tight=False)
                ax2.legend()
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Loop length")
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig1.suptitle('Map of loops')
                fig1.tight_layout()
                return fig1
            else:
                maxG_y = max(loop_listG)
                fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(loop_indG, loop_listG, label=labelF, color=(0,0,1))
                unknown_level=float(maxG_y)*1.1
                first_unknown=True
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1
                    ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax1.autoscale(enable=True, axis='x', tight=False)
                ax1.legend()
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Loop length")
                        
                ax2.plot(0, 0, label=labelR, color=(0,1,0))
                unknown_level= 90
                first_unknown=True
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1
                    ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax2.legend()
                ax2.set_ylim(0, 100)
                ax2.set_xlim(0, xlen)
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Loop length")
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig1.suptitle('Map of loops')
                fig1.tight_layout()
                return fig1

        elif loop_indC is not None:
            maxC_y = max(loop_listC)
            maxC_x = max(loop_indC)
            fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            unknown_level= 90
            first_unknown=True
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1
                ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Loop length")
                        
            ax2.plot(loop_indC, loop_listC, label=labelR, color=(0,1,0))
            unknown_level=float(maxC_y)*1.1
            first_unknown=True
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1
                ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            ax2.autoscale(enable=True, axis='x', tight=False) 
            ax2.legend()
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Loop length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig1.suptitle('Map of loops')
            fig1.tight_layout()
            return fig1
        else:
            fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            unknown_level= 90
            first_unknown=True
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1
                ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Loop length")
                        
            ax2.plot(0, 0, label=labelR, color=(0,1,0))
            unknown_level= 90
            first_unknown=True
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1
                ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            ax2.legend()
            ax2.set_ylim(0, 100)
            ax2.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Loop length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig1.suptitle('Map of loops')
            fig1.tight_layout()
            return fig1


    figure1 = fig1()
    figure1.savefig(os.path.join(save_dir,'plot.png'), dpi=100)
    plt.close(figure1)

    def fig2():

        if loop_list_arrG is not None:
            if loop_list_arrC is not None:
                maxG = np.max(loop_list_arrG)
                maxC = np.max(loop_list_arrC)
                meanG = np.mean(loop_list_arrG)
                meanC = np.mean(loop_list_arrC)
                percentile75G = np.percentile(loop_list_arrG, 75)
                percentile75C = np.percentile(loop_list_arrC, 75)
                percentile25G = np.percentile(loop_list_arrG, 25)
                percentile25C = np.percentile(loop_list_arrC, 25)
                iqrG = percentile75G-percentile25G
                iqrC = percentile75C-percentile25C
                bins1=bins=bin_range(loop_list_arrG, loop_list_arrC)
                maxbin = bins1[-1]
                binGС = int(bins1[1])

                fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.hist([np.clip(loop_list_arrG, bins1[0], bins[-1])], bins=bins1, label=labelF, color=(0,0,1))
                ax1.legend(loc=1)
                plt.xlabel("Loop length")
                ax1.set_xlim(0, maxbin)
                ax1.set_ylabel("Number")
                min_ylim, max_ylim = ax1.set_ylim()
                ax1.axvline(loop_medianG, color='k', linestyle='dashed', linewidth=1)
                ax1.text(loop_medianG*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianG))
                ax1.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanG))
                ax1.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrG))
                ax1.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxG))
                ax1.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
                if maxG >=bins1[-2]:
                    ax1.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                    ax1.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')
                
                ax2.hist([np.clip(loop_list_arrC, bins1[0], bins[-1])], bins=bins1, label=labelR, color=(0,1,0))
                ax2.legend(loc=1)
                plt.xlabel("Loop length")
                ax2.set_xlim(0, maxbin)
                ax2.set_ylabel("Number")
                min_ylim, max_ylim = ax2.set_ylim()
                ax2.axvline(loop_medianC, color='k', linestyle='dashed', linewidth=1)
                ax2.text(loop_medianC*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianC))
                ax2.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanC))
                ax2.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrC))
                ax2.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxC))
                ax2.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
                if maxC >=bins1[-2]:
                    ax2.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                    ax2.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')

                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig2.suptitle('Loop length distribution')
                fig2.tight_layout()
                return fig2
            
            else:
                maxG = np.max(loop_list_arrG)
                meanG = np.mean(loop_list_arrG)
                percentile75G = np.percentile(loop_list_arrG, 75)
                percentile25G = np.percentile(loop_list_arrG, 25)
                iqrG = percentile75G-percentile25G
                bins1=bins=bin_range(loop_list_arrG, loop_list_arrC)
                maxbin = bins1[-1]
                binGС = int(bins1[1])
                fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.hist([np.clip(loop_list_arrG, bins1[0], bins[-1])], bins=bins1, label=labelF, color=(0,0,1))
                ax1.legend(loc=1)
                plt.xlabel("Loop length")
                ax1.set_xlim(0, maxbin)
                ax1.set_ylabel("Number")
                min_ylim, max_ylim = ax1.set_ylim()
                ax1.axvline(loop_medianG, color='k', linestyle='dashed', linewidth=1)
                ax1.text(loop_medianG*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianG))
                ax1.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanG))
                ax1.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrG))
                ax1.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxG))
                ax1.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
                if maxG >=bins1[-2]:
                    ax1.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                    ax1.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')

                ax2.hist(0, bins=1, label=labelR, color=(0,1,0))
                ax2.legend(loc=1)
                plt.xlabel("Loop length")
                ax2.set_xlim(0, maxbin)
                ax2.set_ylabel("Number")
                min_ylim, max_ylim = ax2.set_ylim()
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig2.suptitle('Loop length distribution')
                fig2.tight_layout()
                return fig2
            
        elif loop_list_arrC is not None:
            maxC = np.max(loop_list_arrC)
            meanC = np.mean(loop_list_arrC)
            percentile75C = np.percentile(loop_list_arrC, 75)
            percentile25C = np.percentile(loop_list_arrC, 25)
            iqrC = percentile75C-percentile25C
            bins1=bins=bin_range(loop_list_arrG, loop_list_arrC)
            maxbin = bins1[-1]
            binGС = int(bins1[1])
            fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.hist(0, bins=1, label=labelF, color=(0,0,1))
            ax1.legend(loc=1)
            plt.xlabel("Loop length")
            ax1.set_xlim(0, maxbin)
            ax1.set_ylabel("Number")
            min_ylim, max_ylim = ax1.set_ylim()
            ax2.hist([np.clip(loop_list_arrC, bins1[0], bins[-1])], bins=bins1, label=labelR, color=(0,1,0))
            ax2.legend(loc=1)
            plt.xlabel("Loop length")
            ax2.set_xlim(0, maxbin)
            ax2.set_ylabel("Number")
            min_ylim, max_ylim = ax2.set_ylim()
            ax2.axvline(loop_medianC, color='k', linestyle='dashed', linewidth=1)
            ax2.text(loop_medianC*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianC))
            ax2.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanC))
            ax2.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrC))
            ax2.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxC))
            ax2.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
            if maxC >=bins1[-2]:
                ax2.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                ax2.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')
            
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig2.suptitle('Loop length distribution')
            fig2.tight_layout()
            return fig2
        
        else:
            maxbin = 2000
        
            fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.hist(0, bins=1, label=labelF, color=(0,0,1))
            ax1.legend(loc=1)
            plt.xlabel("Loop length")
            ax1.set_xlim(0, maxbin)
            ax1.set_ylabel("Number")
            min_ylim, max_ylim = ax1.set_ylim()

            ax2.hist(0, bins=1, label=labelR, color=(0,1,0))
            ax2.legend(loc=1)
            plt.xlabel("Loop length")
            ax2.set_xlim(0, maxbin)
            ax2.set_ylabel("Number")
            min_ylim, max_ylim = ax2.set_ylim()
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig2.suptitle('Loop length distribution')
            fig2.tight_layout()
            return fig2


    figure2 = fig2()
    figure2.savefig(os.path.join(save_dir, 'plot_hist.png'), dpi=100)
    plt.close(figure2)

    # Derive session folder name
    session_folder = os.path.basename(save_dir)
    # Create URLs compounds that Flask can serve
    plot1_url = f"{session_folder}/plot.png"
    plot2_url = f"{session_folder}/plot_hist.png"
    return plot1_url, plot2_url  

def repeats_dataframeGC(id_, ind_start, ind_end):

    tel_df = pd.DataFrame({   
        "Sequence name": [id_]*len(ind_start),
        "Tel Start": ind_start,
        "Tel End" : ind_end})
    
    return tel_df

def csvGC_repeats_open(patternG, patternC, save_dir, selected_option):
        
    csv_path_repeatsG = os.path.join(save_dir, f'repeatsF.csv')
    csv_file_repeatsG = open(csv_path_repeatsG, 'w', newline='')
    headersG_repeats = ['Sequence name', 'Start', 'End', f'Pattern searched {patternG}']
    writer = csv.writer(csv_file_repeatsG, delimiter=',')
    writer.writerow(headersG_repeats)
    
    csv_path_repeatsC = os.path.join(save_dir, f'repeatsR.csv')
    csv_file_repeatsC = open(csv_path_repeatsC, 'w', newline='')
    headersC_repeats = ['Sequence name', 'Start', 'End', f'Pattern searched {patternC}']
    writer = csv.writer(csv_file_repeatsC, delimiter=',')
    writer.writerow(headersC_repeats)

    return csv_path_repeatsG, csv_path_repeatsC

def save_repeats_csvGC(repeatsG, repeatsC, csv_file_repeatsG, csv_file_repeatsC):
    repeatsG.to_csv(csv_file_repeatsG, mode='a', sep=',', index=False, header=False, date_format=object)
    repeatsC.to_csv(csv_file_repeatsC, mode='a', sep=',', index=False, header=False, date_format=object)

    return repeatsG, repeatsC

def csvGC_clusters_open(patternG, patternC, clusterscore, input_looplength, input_ssr, input_repeats,save_dir, selected_option):

    csv_path_clustersG = os.path.join(save_dir, f'clustersF_score_{clusterscore}_loop_{input_looplength}_ssr_{input_ssr}_nrepeats_{input_repeats}.csv')
    csv_file_clustersG = open(csv_path_clustersG, 'w', newline='')
    headersG = ['Sequence name', 'Cluster Start', 'Cluster End', 'Cluster length', 'Number of repeats', 'Total Loops length', 'Score', 'Input loop', 'SSR', f'Pattern searched {patternG}']
    writer = csv.writer(csv_file_clustersG, delimiter=',')
    writer.writerow(headersG)

    csv_path_clustersC = os.path.join(save_dir, f'clustersR_score_{clusterscore}_loop_{input_looplength}_ssr_{input_ssr}_nrepeats_{input_repeats}.csv')
    csv_file_clustersC = open(csv_path_clustersC, 'w', newline='')
    headersC = ['Sequence name', 'Cluster Start', 'Cluster End', 'Cluster length', 'Number of repeats', 'Total Loops length', 'Score', 'Input loop', 'SSR', f'Pattern searched {patternC}']
    writer = csv.writer(csv_file_clustersC, delimiter=',')
    writer.writerow(headersC)
    return csv_path_clustersG, csv_path_clustersC

    
def save_clusters_csvGC(clusters_dfG, clusters_dfC, csv_file_clustersG, csv_file_clustersC):
    clusters_dfG.to_csv(csv_file_clustersG, mode='a', sep=',', index=False, header=False, date_format=object)
    clusters_dfC.to_csv(csv_file_clustersC, mode='a', sep=',', index=False, header=False, date_format=object)
    return csv_file_clustersG, csv_file_clustersC


def image_loops_multifasta(myseq, loop_indG, loop_listG, loop_indC, loop_listC, description, loops_map, selected_option):
    """
    Loops map for multifasta

    """
    if selected_option == 'TTAGGG':
        labelF = 'LoopG'
        labelR = 'LoopC'
    
    elif selected_option == 'FuzzyTel':
        labelF = 'LoopG'
        labelR = 'LoopC'

    elif selected_option == 'Custom pattern':
        labelF = 'LoopF'
        labelR = 'LoopR'
    
    set_fontsize = 8
    xlen = len(myseq)

    # degenerate_patterns='[NRYBDKMHVSW]+'
    degenerate_patterns = '[^ATGC]+'

    unknown_pattern=re.compile(degenerate_patterns,re.IGNORECASE)
    unknown_ind=[( x.start(), x.end() ) for x in re.finditer(unknown_pattern,myseq)]

    def fig1():
        if loop_indG is not None:
            if loop_indC is not None:
                maxG_y = max(loop_listG)
                maxC_y = max(loop_listC)
                fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(loop_indG, loop_listG, label=labelF, color=(0,0,1))
                unknown_level=float(maxG_y)*1.1
                first_unknown=True
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1  
                    ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax1.legend()
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Loop length")
                        
                ax2.plot(loop_indC, loop_listC, label=labelR, color=(0,1,0))
                unknown_level=float(maxC_y)*1.1
                first_unknown=True
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1  
                    ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax2.legend()
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Loop length")
                ax1.autoscale(enable=True, axis='x', tight=False)
                ax2.autoscale(enable=True, axis='x', tight=False)
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig1.suptitle('Map of loops')
                return fig1
            else:
                maxG_y = max(loop_listG)
                fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(loop_indG, loop_listG, label=labelF, color=(0,0,1))
                unknown_level=float(maxG_y)*1.1
                first_unknown=True
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1  
                    ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax1.legend()
                ax1.autoscale(enable=True, axis='x', tight=False)
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Loop length")
                        
                ax2.plot(0, 0, label=labelR, color=(0,1,0))
                first_unknown=True
                unknown_level=90
                for unknown_start,unknown_end in unknown_ind:
                    if abs(unknown_end - unknown_start) < 2: 
                        unknown_end = unknown_end + 1  
                    ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                    first_unknown=False
                ax2.legend()
                ax2.set_ylim(0, 100)
                ax2.set_xlim(0, xlen)
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Loop length")
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig1.suptitle('Map of loops')
                return fig1

        elif loop_indC is not None: 
            maxC_y = max(loop_listC)  
            fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            first_unknown=True
            unknown_level=90
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1  
                ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Loop length")
                        
            ax2.plot(loop_indC, loop_listC, label=labelR, color=(0,1,0))
            unknown_level=float(maxC_y)*1.1
            first_unknown=True
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1  
                ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            ax2.legend()
            ax2.autoscale(enable=True, axis='x', tight=False)
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Loop length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig1.suptitle('Map of loops')
            return fig1
        else:
            fig1, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            first_unknown=True
            unknown_level=90
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1  
                ax1.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Loop length")
                        
            ax2.plot(0, 0, label=labelR, color=(0,1,0))
            unknown_level=90
            first_unknown=True
            for unknown_start,unknown_end in unknown_ind:
                if abs(unknown_end - unknown_start) < 2: 
                    unknown_end = unknown_end + 1  
                ax2.plot([unknown_start,unknown_end],[unknown_level,unknown_level],label='Redundant' if first_unknown else None,color=(1,0,0))
                first_unknown=False
            ax2.legend()
            ax2.set_ylim(0, 100)
            ax2.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Loop length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig1.suptitle('Map of loops')
            return fig1

    
    figure1 = fig1()
    loops_map.append(figure1)
    plt.close(figure1)
    return loops_map

def image_hist_loops_multifasta(loop_list_arrG, loop_list_arrC, description, loops_hist, loop_medianG, loop_medianC, selected_option):
    """
    Hist map for multifasta

    """
    if selected_option == 'TTAGGG':
        labelF = 'LoopG'
        labelR = 'LoopC'
    
    elif selected_option == 'FuzzyTel':
        labelF = 'LoopG'
        labelR = 'LoopC'

    elif selected_option == 'Custom pattern':
        labelF = 'LoopF'
        labelR = 'LoopR'
    
    set_fontsize = 8

    def fig2():

        if loop_list_arrG is not None:
            if loop_list_arrC is not None:
                maxG = np.max(loop_list_arrG)
                maxC = np.max(loop_list_arrC)
                meanG = np.mean(loop_list_arrG)
                meanC = np.mean(loop_list_arrC)
                percentile75G = np.percentile(loop_list_arrG, 75)
                percentile75C = np.percentile(loop_list_arrC, 75)
                percentile25G = np.percentile(loop_list_arrG, 25)
                percentile25C = np.percentile(loop_list_arrC, 25)
                iqrG = percentile75G-percentile25G
                iqrC = percentile75C-percentile25C
                bins1=bins=bin_range(loop_list_arrG, loop_list_arrC)
                maxbin = bins1[-1]
                binGС = int(bins1[1])

                fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.hist([np.clip(loop_list_arrG, bins1[0], bins[-1])], bins=bins1, label=labelF, color=(0,0,1))
                ax1.legend(loc=1)
                plt.xlabel("Loop length")
                ax1.set_xlim(0, maxbin)
                ax1.set_ylabel("Number")
                min_ylim, max_ylim = ax1.set_ylim()
                ax1.axvline(loop_medianG, color='k', linestyle='dashed', linewidth=1)
                ax1.text(loop_medianG*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianG))
                ax1.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanG))
                ax1.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrG))
                ax1.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxG))
                ax1.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
                if maxG >=bins1[-2]:
                    ax1.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                    ax1.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')
                
                ax2.hist([np.clip(loop_list_arrC, bins1[0], bins[-1])], bins=bins1, label=labelR, color=(0,1,0))
                ax2.legend(loc=1)
                plt.xlabel("Loop length")
                ax2.set_xlim(0, maxbin)
                ax2.set_ylabel("Number")
                min_ylim, max_ylim = ax2.set_ylim()
                ax2.axvline(loop_medianC, color='k', linestyle='dashed', linewidth=1)
                ax2.text(loop_medianC*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianC))
                ax2.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanC))
                ax2.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrC))
                ax2.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxC))
                ax2.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
                if maxC >=bins1[-2]:
                    ax2.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                    ax2.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig2.suptitle('Loop length distribution')
                return fig2

                
            else:
                maxG = np.max(loop_list_arrG)
                meanG = np.mean(loop_list_arrG)
                percentile75G = np.percentile(loop_list_arrG, 75)
                percentile25G = np.percentile(loop_list_arrG, 25)
                iqrG = percentile75G-percentile25G
                bins1=bins=bin_range(loop_list_arrG, loop_list_arrC)
                maxbin = bins1[-1]
                binGС = int(bins1[1])

                fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.hist([np.clip(loop_list_arrG, bins1[0], bins[-1])], bins=bins1, label=labelF, color=(0,0,1))
                ax1.legend(loc=1)
                plt.xlabel("Loop length")
                ax1.set_xlim(0, maxbin)
                ax1.set_ylabel("Number")
                min_ylim, max_ylim = ax1.set_ylim()
                ax1.axvline(loop_medianG, color='k', linestyle='dashed', linewidth=1)
                ax1.text(loop_medianG*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianG))
                ax1.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanG))
                ax1.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrG))
                ax1.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxG))
                ax1.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
                if maxG >=bins1[-2]:
                    ax1.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                    ax1.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')
                
                ax2.hist(0, bins=1, label=labelR, color=(0,1,0))
                ax2.legend()
                plt.xlabel("Loop length")
                ax2.set_xlim(0, maxbin)
                ax2.set_ylabel("Number")
                min_ylim, max_ylim = ax2.set_ylim()
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig2.suptitle('Loop length distribution')
                return fig2
            
        elif loop_list_arrC is not None:
            maxC = np.max(loop_list_arrC)
            meanC = np.mean(loop_list_arrC)
            percentile75C = np.percentile(loop_list_arrC, 75)
            percentile25C = np.percentile(loop_list_arrC, 25)
            iqrC = percentile75C-percentile25C
            bins1=bins=bin_range(loop_list_arrG, loop_list_arrC)
            maxbin = bins1[-1]
            binGС = int(bins1[1])

            fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.hist(0, bins=1, label=labelF, color=(0,0,1))
            ax1.legend()
            plt.xlabel("Loop length")
            ax1.set_xlim(0, maxbin)
            ax1.set_ylabel("Number")
            min_ylim, max_ylim = ax1.set_ylim()
            
            ax2.hist([np.clip(loop_list_arrC, bins1[0], bins[-1])], bins=bins1, label=labelR, color=(0,1,0))
            ax2.legend(loc=1)
            plt.xlabel("Loop length")
            ax2.set_xlim(0, maxbin)
            ax2.set_ylabel("Number")
            min_ylim, max_ylim = ax2.set_ylim()
            ax2.axvline(loop_medianC, color='k', linestyle='dashed', linewidth=1)
            ax2.text(loop_medianC*1.1, max_ylim*0.8, 'Median:\n {:.2f}'.format(loop_medianC))
            ax2.text(maxbin*0.75, max_ylim*0.7, 'Mean: {:.2f}'.format(meanC))
            ax2.text(maxbin*0.75, max_ylim*0.6, 'IQR: {:.2f}'.format(iqrC))
            ax2.text(maxbin*0.75, max_ylim*0.5, 'Max: {:d}'.format(maxC))
            ax2.text(maxbin*0.75, max_ylim*0.4, 'Bin: {:d}'.format(binGС))
            if maxC >=bins1[-2]:
                ax2.axvline(bins1[-2], color='r', linestyle='dashed', linewidth=1)
                ax2.text(maxbin*0.75, max_ylim*0.3, 'OL: \u2265 {:d}'.format(math.floor(bins1[-2])), color='r')
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig2.suptitle('Loop length distribution')
            return fig2
        
        else:
            fig2, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.hist(0, bins=1, label=labelF, color=(0,0,1))
            ax1.legend()
            plt.xlabel("Loop length")
            ax1.set_xlim(0, 2000)
            ax1.set_ylabel("Number")
            min_ylim, max_ylim = ax1.set_ylim()


            ax2.hist(0, bins=1, label=labelR, color=(0,1,0))
            ax2.legend()
            plt.xlabel("Loop length")
            ax2.set_xlim(0, 2000)
            ax2.set_ylabel("Number")
            min_ylim, max_ylim = ax2.set_ylim()
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig2.suptitle('Loop length distribution')
            return fig2


    figure2 = fig2()
    loops_hist.append(figure2)
    plt.close(figure2)
    return loops_hist

def image_clusters_multifasta(myseq, clusterG_list, clusterG_ind, clusterC_list, clusterC_ind, description, clusters_map, selected_option):
    """
    Cluster map for multifasta

    """
    if selected_option == 'TTAGGG':
        labelF = 'ClusterG'
        labelR = 'ClusterC'
    
    elif selected_option == 'FuzzyTel':
        labelF = 'ClusterG'
        labelR = 'ClusterC'

    elif selected_option == 'Custom pattern':
        labelF = 'ClusterF'
        labelR = 'ClusterR'
    
    set_fontsize = 8

    xlen = len(myseq)
    
    def fig3():
        if sum(clusterG_list) > 0:
            if sum(clusterC_list) > 0:
                fig3, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(clusterG_ind, clusterG_list, label=labelF, color=(0,0,1))
                ax1.legend()
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Cluster length")
                            
                ax2.plot(clusterC_ind, clusterC_list, label=labelR, color=(0,1,0))
                ax2.legend()
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Cluster length")
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig3.suptitle('Map of Repeat clusters')
                return fig3
            else:
                fig3, (ax1, ax2) = plt.subplots(2,1, sharex=True)
                ax1.set_title(description, fontsize=set_fontsize)
                ax1.plot(clusterG_ind, clusterG_list, label=labelF, color=(0,0,1))
                ax1.legend()
                plt.xlabel("Sequence position")
                ax1.set_ylabel("Cluster length")
                            
                ax2.plot(0, 0, label=labelR, color=(0,1,0))
                ax2.legend()
                ax2.set_ylim(0, 100)
                ax2.set_xlim(0, xlen)
                plt.xlabel("Sequence position")
                ax2.set_ylabel("Cluster length")
                ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
                ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
                fig3.suptitle('Map of Repeat clusters')
                return fig3

        elif sum(clusterC_list) > 0:    
            fig3, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Cluster length")
                            
            ax2.plot(clusterC_ind, clusterC_list, label=labelR, color=(0,1,0))
            ax2.legend()
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Cluster length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig3.suptitle('Map of Repeat clusters')
            return fig3
        else:
            fig3, (ax1, ax2) = plt.subplots(2,1, sharex=True)
            ax1.set_title(description, fontsize=set_fontsize)
            ax1.plot(0, 0, label=labelF, color=(0,0,1))
            ax1.legend()
            ax1.set_ylim(0, 100)
            ax1.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax1.set_ylabel("Cluster length")
                            
            ax2.plot(0, 0, label=labelR, color=(0,1,0))
            ax2.legend()
            ax2.set_ylim(0, 100)
            ax2.set_xlim(0, xlen)
            plt.xlabel("Sequence position")
            ax2.set_ylabel("Cluster length")
            ax1.xaxis.set_minor_locator(AutoMinorLocator(10))
            ax2.xaxis.set_minor_locator(AutoMinorLocator(10)) 
            ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax2.yaxis.set_minor_locator(AutoMinorLocator(4)) 
            fig3.suptitle('Map of Repeat clusters')
            return fig3

    figure3 = fig3()
    clusters_map.append(figure3)
    plt.close(figure3)
    return clusters_map


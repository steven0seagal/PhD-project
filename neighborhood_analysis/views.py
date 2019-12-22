from django.shortcuts import render, redirect
from django.http import HttpResponse
from accounts.models import CompleteQueue
from .models import NeighAnalyzDatabase
from scripts.create_cron_task import CreateCronTaskNeighborhood
import random
import string
from scripts.PfamValidate import PfamValidator
from scripts.pfam_input import retype_domain, check_if_domain_can_be_reachable

def neighana(request):
    return render(request, 'tools/input/neighborhood_analyzer.html')

def count(request):
#    fulltext = request.GET['fulltext']
    pfam_domain = request.GET['pfam_domain']
    range_search = request.GET['range_search']
    cut_off = request.GET['cut_off']
    database_taxa = request.GET['database_taxa']
    strand_select = request.GET['strand_select']
    out_name = request.GET['out_name']
    test_correction = request.GET['test_correction'] 
    #out_time = '12:34:56.789012'
    user_id  = request.user.id


    """
    # Beauty for dashboard
    if database_taxa == 'alldb':
        smart_database_taxa = 'All database'
    else:
        smart_database_taxa = database_taxa.capitalize() + ' sp.'
    if strand_select == 'same':
        smart_strand_select = 'Same as searched domain '
    elif strand_select == 'oposite':
        smart_strand_select = 'Oposite as searched domain'
    else:
        smart_strand_select = 'Both strands'
    if test_correction == 'fdr_bh':
        smart_test_correction = 'Benjamini-Hochberg Procedure'
    elif test_correction == 'bonferroni':
        smart_test_correction = 'Bonferroni Correction'
    else:
        smart_test_correction = 'None' 
    """
    pfam = PfamValidator(pfam_domain)
    pfam_value = pfam.Validate()
    domain_for_database, domain_to_search = retype_domain(pfam_domain)

    if pfam_value is True:
        if check_if_domain_can_be_reachable(database_taxa,domain_to_search) is True:
            letters = string.ascii_lowercase
            end_end =  ''.join(random.choice(letters) for i in range(15))
            link_down = '/media/results/' + end_end +'.txt'
            # Local
            # ready_script ='python3 /mnt/d/45.76.38.24/final_site-project/scripts/execute_order_66.py {} {} {} ' \
            #              '{} {} {} {} '.format(domain_to_search ,range_search, database_taxa, cut_off, test_correction,
            #                                 strand_select, link_down)
            # Vultr
            ready_script = '/home/djangoadmin/final_site_venv/bin/python3 /home/djangoadmin/final_site-project/scripts/execute_order_66.py {} {} {} {} {} {} {}'.format(domain_to_search ,range_search, database_taxa, cut_off, test_correction,strand_select, link_down)
                            

            job = CompleteQueue(user_id = user_id,tool = 'NA', status = 'Queue', analysis_name = out_name,
                                script = ready_script, file = link_down )
            job.save()

            return redirect('dashboard')
        else:
            return render(request, 'tools/error/heavy_calculation.html')
    else:
        return render(request,'tools/error/wrong_pfam.html')

"""
    a = NeighbourhoodAnalyzer(pfam_domain,  range_search, database_taxa)
    complete_data, message_down,message_for_output,message_for_additional_data = a.GO()
    
    context ={
        'complete_data' :  complete_data,
        'message_down' : message_down,
        'message_for_output' : message_for_output,
        'message_for_additional_data' : message_for_additional_data
            }
    #link_ready = results_link+output_file+'.txt'
    


    return render(request,'tools/output/neighborhood_result.html',context )
"""

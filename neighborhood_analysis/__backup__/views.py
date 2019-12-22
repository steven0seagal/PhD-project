from django.shortcuts import render, redirect
from django.http import HttpResponse
from scripts.analiza_otoczenia_class_v2_0_linux import NeighbourhoodAnalyzer
from .models import NeighAnalyzDatabase
from scripts.create_cron_task import CreateCronTaskNeighborhood
from datetime import datetime, timedelta
from scripts.PfamValidate import PfamValidator
from scripts.pfam_input import retype_domain

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
    #task = CreateCronTask(pfam_domain, int(range_search), database_taxa)
    #out_time = task.Task()
    out_time = '12:34:56.789012'
    #link_down = '/home/site/results/'+str(out_time)+'.txt'
    #send info to database 
    user_id  = request.user.id
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

    pfam = PfamValidator(pfam_domain)
    pfam_value = pfam.Validate()
    domain_for_database, domain_to_search = retype_domain(pfam_domain)
    if pfam_value == True:

        create_task = CreateCronTaskNeighborhood(domain_to_search, range_search, database_taxa, cut_off,test_correction,strand_select, out_name)
        pre_end = create_task.Task()
        end_end = pre_end.replace(' ','_')
        link_down = '/media/' + end_end +'.txt'
        query = NeighAnalyzDatabase(pfam_name = domain_for_database, neigh_size = range_search, cut_off = cut_off, tax = smart_database_taxa, strand = smart_strand_select, user_id = user_id, out_name = out_name, out_time = end_end ,test_correction = smart_test_correction,link = link_down )
        query.save()

        return redirect('dashboard')
    elif pfam_value == False:
         
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

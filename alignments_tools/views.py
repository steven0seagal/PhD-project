from django.shortcuts import render,redirect
from django.http import HttpResponse
from scripts.colapsing_script import ColapsingSeq
from django.core.files.storage import FileSystemStorage
from scripts.streching3alignments import ThreeAlignmentColapser
from .models import ColapserDatabase, StretcherDatabase
from scripts.create_cron_task import CreateCronTaskStreching
from accounts.models import CompleteQueue
import string,random

def match3align(request):
    return render(request, 'tools/input/match_3_alignment.html')


def squezealign(request):
    return render(request, 'tools/input/colapsing_alignment.html')




def makeitsimple(request):
    context = {}
    #WINDOWS
    #linki = 'E:/project_site/final_site-project'
    #LINUX
    linki = '/home/djangoadmin/final_site-project'
    results_link = '/media/'
    
    if request.user.is_authenticated:
        
        if request.method =='POST':
            out_name = request.POST['out_name']
            uploaded_file = request.FILES['document']
            fs = FileSystemStorage()
            name= fs.save(uploaded_file.name, uploaded_file)
            context['url'] = fs.url(name)
            link_do_pliku = context['url']


            a = ColapsingSeq(multifasta_sequences = linki + link_do_pliku)
            output = a.Go()
            guide_to_file = a.CreateLink()
            
            
            link_ready = results_link+guide_to_file  
            user_id  = request.user.id
            query = ColapserDatabase(insert_file = name, user_id = user_id, out_time = 'NOW', link = link_ready,anal_name = out_name)
            query.save()

            return render(request, 'tools/output/colapsed.html', {'output':output, 'link_ready':link_ready})

    else:
         if request.method =='POST':
            out_name = request.POST['out_name']
            uploaded_file = request.FILES['document']
            fs = FileSystemStorage()
            name= fs.save(uploaded_file.name, uploaded_file)
            context['url'] = fs.url(name)
            link_do_pliku = context['url']


            a = ColapsingSeq(multifasta_sequences = linki + link_do_pliku)
            output = a.Go()
            guide_to_file = a.CreateLink()       
            
             
            link_ready = results_link+guide_to_file  
            user_id  = request.user.id
            query = ColapserDatabase(insert_file = name, user_id = 2, out_time = 'NEVER', link = link_ready,anal_name = out_name)
            query.save()



            return render(request, 'tools/output/colapsed.html', {'output':output, 'link_ready':link_ready})




def strechitout(request):
    #Skoro to jest slownik to mo¿na to za³atwiæ jednym a nie trzema
    context_main = {}
    context_subordinate1 = {}
    context_subordinate2 = {}
    linki = '/home/djangoadmin/final_site-project/media/'
    
    results_link = '/media/results/'
    if request.method =='POST':
        out_name = request.POST['out_name']
        uploaded_file_main = request.FILES['main']
        uploaded_file_subordinate1 = request.FILES['subordinate1']
        uploaded_file_subordinate2 = request.FILES['subordinate2']
        fs = FileSystemStorage()
        name_main = fs.save(uploaded_file_main.name, uploaded_file_main)
        name_subordinate1 = fs.save(uploaded_file_subordinate1.name, uploaded_file_subordinate1)
        name_subordinate2 = fs.save(uploaded_file_subordinate2.name, uploaded_file_subordinate2)
        context_main['url'] = fs.url(name_main)
        context_subordinate1['url'] = fs.url(name_subordinate1)
        context_subordinate2['url'] = fs.url(name_subordinate2)
        link_do_main = context_main['url']
        link_do_subordinate1 = context_subordinate1['url']
        link_do_subordinate2 = context_subordinate2['url']
        user_id = request.user.id

        letters = string.ascii_lowercase
        end_end = ''.join(random.choice(letters) for i in range(15))
        link_down ='/media/results/' + end_end + '_streched.txt'
        ################################################################################################################
        # Local
        ready_script = 'python3 /mnt/d/45.76.38.24/final_site-project/scripts/execute_order_69.py {} {} ' \
                       '{} {}'.format(name_main, name_subordinate1,name_subordinate2, link_down)
        # Vultr
        # ready_script = """ '/home/djangoadmin/final_site_venv/bin/python3
        #                   /home/djangoadmin/final_site-project/scripts/execute_order_69.py {} {} ' \
        #                   {} {}'.format(name_main, name_subordinate1,name_subordinate2, link_down)"""

        job = CompleteQueue(user_id=user_id, tool='M3A', status='Queue', analysis_name=out_name,
                            script=ready_script, file=link_down)
        job.save()

        ################################################################################################################
        """
        create_task = CreateCronTaskStreching(name_main, name_subordinate1, name_subordinate2, out_name)
        pre_end = create_task.Task()
        end_end = pre_end.replace(' ','_')
        link_down = '/media/' + end_end +'.txt'
        #a = ThreeAlignmentColapser(reference_multifasta = name_main , first_multifasta = name_subordinate1, second_multifasta = name_subordinate2)
        #a.Go()
        #guide_to_file = a.CreateDownloadableLink()
        #link_ready = results_link+guide_to_file
        query = StretcherDatabase(insert_file_main_main = name_main , insert_file_main_subordinate_one = name_subordinate1, insert_file_main_subordinate_two= name_subordinate2, user_id = user_id, out_time = end_end, link = link_down, anal_name = out_name)
        query.save()
        """
    return redirect('dashboard')
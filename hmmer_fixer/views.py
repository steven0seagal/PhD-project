from django.shortcuts import render, redirect
from django.http import HttpResponse
from .models import HmmerFixerDatabase
from django.core.files.storage import FileSystemStorage
from scripts.HMMER_seq_fix_BB import HMMER


# Create your views here.

def hmmer_start(request):
    title = " HMMER fixer v.0.1 "
    return render(request, 'tools/input/hmmer_fixer_start.html')


def hmmer_analysis(request):
    context_main = {}
    context_skipped = {}

    # WINDOWS
    # linki = 'E:/project_site/final_site-project'
    # LINUX
    linki = '/home/djangoadmin/final_site-project'

    results_link = '/media/'

    if request.user.is_authenticated:

        if request.method == 'POST':
            fs = FileSystemStorage()
            out_name = request.POST['out_name']

            uploaded_file_main = request.FILES['master']
            name_main = fs.save(uploaded_file_main.name, uploaded_file_main)
            context_main['url'] = fs.url(name_main)
            link_do_main = context_main['url']
            try:
                uploaded_file_skipped = request.FILES['skipdaseks']
                name_skipped = fs.save(uploaded_file_skipped.name, uploaded_file_skipped)
                context_skipped['url'] = fs.url(name_skipped)
                link_do_skipped = context_skipped['url']
                a = HMMER(main_seqs=linki + link_do_main, skip_seqs=linki + link_do_skipped)
            except:
                name_skipped = "None"
                a = HMMER(main_seqs=linki + link_do_main, skip_seqs=None)

            a.Go()
            guide_to_file = a.result_name_maker()
            guide_to_file = guide_to_file.split('/')
            link_ready = results_link + guide_to_file[-1]
            user_id = request.user.id
            query = HmmerFixerDatabase(main_file=name_main, skipped_seqs=name_skipped, user_id=user_id,
                                       anal_name=out_name, link=link_ready)
            query.save()
            data = "You ve dane everything correct"
            content = {
                'data': data,
                'link_ready': link_ready
            }
            return render(request, 'tools/output/hmmer_fixed.html', content)

    else:
        if request.method == "POST":
            fs = FileSystemStorage()
            out_name = request.POST['out_name']

            uploaded_file_main = request.FILES['master']
            name_main = fs.save(uploaded_file_main.name, uploaded_file_main)
            context_main['url'] = fs.url(name_main)
            link_do_main = context_main['url']
            try:
                uploaded_file_skipped = request.FILES['skipdaseks']
                name_skipped = fs.save(uploaded_file_skipped.name, uploaded_file_skipped)
                context_skipped['url'] = fs.url(name_skipped)
                link_do_skipped = context_skipped['url']
                a = HMMER(main_seqs=linki + link_do_main, skip_seqs=linki + link_do_skipped)
            except:
                name_skipped = "None"
                a = HMMER(main_seqs=linki + link_do_main, skip_seqs=None)
            a.Go()
            guide_to_file = a.result_name_maker()
            guide_to_file = guide_to_file.split('/')
            link_ready = results_link + guide_to_file[-1]



            query = HmmerFixerDatabase(main_file=name_main, skipped_seqs=name_skipped, user_id=5, anal_name=out_name,
                                       link='testowy_down')
            query.save()

            data = "You ve dane everything correct"

            content = {
                'data': data,
                'link_ready': link_ready,
            }
            return render(request, 'tools/output/hmmer_fixed.html', content)
from django.shortcuts import render, redirect
from django.contrib import messages, auth
from django.contrib.auth.models import User
from neighborhood_analysis.models import NeighAnalyzDatabase
from alignments_tools.models import ColapserDatabase, StretcherDatabase
from hmmer_fixer.models import HmmerFixerDatabase
from datetime import datetime
from scripts import check_link_existance
from .models import PeekUserData
from .models import CompleteQueue
from scripts.data_from_script import feedMe
from scripts.check_link_existance import new_file_checker
import os
def register(request):
    if request.method == "POST":
        #Get form values
        first_name = request.POST['first_name']
        last_name = request.POST['last_name']
        username = request.POST['username']
        email = request.POST['email']
        password = request.POST['password']
        password2 = request.POST['password2']

        #Check if password match
        if password == password2:
            #check username
            if User.objects.filter(username=username).exists():
                messages.error(request, 'That username is taken')
                return redirect('register')
            else:
                if User.objects.filter(email=email).exists():
                    messages.error(request, 'That email is being used')
                    return redirect('register')
                else:
                    user = User.objects.create_user(username=username,email=email, password=password, first_name=first_name, last_name=last_name)
                    query = PeekUserData(username=username,email=email, password=password, first_name=first_name, last_name=last_name)
                    #Login after register 
                    # auth.login(request, user)
                    # messages.success(request, 'You are now registered')
                    # return redirect('index')
                    user.save()
                    query.save()
                    messages.success(request, 'You are now registered and can log in ')
                    return redirect('login')
        else:
            messages.error(request, 'Passwords do not match')
            return redirect('register')
    else:
        return render(request, 'user_panel/register.html')

def login(request):
    if request.method == "POST":
        username=request.POST['username']
        password=request.POST['password']

        user = auth.authenticate(username=username, password=password)
       
        if user is not None:
            auth.login(request, user)
            messages.success(request, 'You are now logged in')
            return redirect('dashboard')
        else:
            messages.error(request, 'Invalid credentials')
            return redirect('login')
    else:
        return render(request, 'user_panel/login.html')

def logout(request):
    if request.method == "POST":
        auth.logout(request)
        messages.success(request, 'You are now logout')

        return redirect('login')

def dashboard(request):
    current_time =  str(datetime.now())


    colapser_jobs = ColapserDatabase.objects.order_by('-id').filter(user_id = request.user.id)
    stretcher_jobs = StretcherDatabase.objects.order_by('-id').filter(user_id = request.user.id)
    hmmer_jobs = HmmerFixerDatabase.objects.order_by('-id').filter(user_id = request.user.id)

    all_jobs = CompleteQueue.objects.order_by('-id').filter(user_id = request.user.id)
    neigh_example_jobs = feedMe(all_jobs, 'NA')
    stretch_example_jobs = feedMe(all_jobs,'M3A')




    checked_data_colapser = check_link_existance.link_ready(colapser_jobs)
    checked_strecher_jobs = check_link_existance.link_ready(stretcher_jobs)
    checked_hmmer_jobs = check_link_existance.link_ready(hmmer_jobs)


    context = {
	    'current_time' : current_time,
        'stretch_example_jobs' :stretch_example_jobs,
        'checked_data_colapser' : checked_data_colapser,
        'checked_hmmer_jobs' : checked_hmmer_jobs,
        'neigh_example_jobs':neigh_example_jobs,
    }

#    context = {
#        'neigh_jobs' : neigh_jobs,
#        'colapser_jobs' : colapser_jobs,
#        'stretcher_jobs' : stretcher_jobs,
#	    'current_time' : current_time,
#        'hmmer_jobs' : hmmer_jobs,
#    }
    return render(request, 'user_panel/dashboard.html',context)
from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse
from.models import Tools
from django.core.paginator import Paginator

def mainsite(request):
    tools_list = Tools.objects.all()
    paginator = Paginator(tools_list, 6)
    page = request.GET.get('page')
    tools = paginator.get_page(page)
    return render(request, 'mainsite/index.html',{'tools':tools})


from django.shortcuts import render, get_object_or_404
from .models import NewsRoom
from django.core.paginator import Paginator

# Create your views here.
def allnews(request):
    news_list = NewsRoom.objects.all()
    paginator = Paginator(news_list, 5)
    page = request.GET.get('page')
    news = paginator.get_page(page)
    return render(request, 'news/newest.html', {'news':news})

def detail(request,news_id):
    detailnews = get_object_or_404(NewsRoom, pk= news_id)
    return render(request, 'news/detail.html',{'news': detailnews})
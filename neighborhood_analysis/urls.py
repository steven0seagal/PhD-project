from django.urls import path
from . import views

urlpatterns = [
    path('', views.neighana_menu, name = 'neighana_menu'),
    path('count', views.count, name = 'count'),
    path('pfam_domain',views.neighana, name = 'neighana'),
    path('gene', views.geneana, name = 'geneana')
]
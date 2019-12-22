from django.urls import path
from . import views

urlpatterns = [
    path('',views.hmmer_start, name = 'hmmer_start'),
    path('hmmer_analysis', views.hmmer_analysis,name = 'hmmer_analysis' ),
]
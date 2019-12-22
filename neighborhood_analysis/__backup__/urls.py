from django.urls import path
from . import views

urlpatterns = [
    path('', views.neighana, name = 'neighana'),
    path('count', views.count, name = 'count'),

]
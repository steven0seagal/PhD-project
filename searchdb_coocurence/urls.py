from django.urls import path
from . import views

urlpatterns = [
    path('', views.searchdbcooc, name = 'searchdbcooc'),
    path('result', views.searchresult, name = 'searchresult'),
]
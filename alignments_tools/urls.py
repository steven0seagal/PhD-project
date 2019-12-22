from django.urls import path
from . import views
urlpatterns = [
    path('match3align', views.match3align, name = 'match3align'), 
    path('squezealign', views.squezealign, name = 'squezealign'), 
    path('makeitsimple', views.makeitsimple, name='makeitsimple'),
    path('strechitout', views.strechitout, name='strechitout'),
    
    
]
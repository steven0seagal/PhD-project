from django.urls import path

from . import views


urlpatterns = [
    path('',views.allnews, name = 'allnews'),
    path('<int:news_id>/', views.detail, name = 'detail')
]
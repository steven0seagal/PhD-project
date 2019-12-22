from django.contrib import admin

# Register your models here.
from .models import ColapserDatabase, StretcherDatabase

class ColapserDatabaseAdmin(admin.ModelAdmin):
    list_display = ('anal_name','insert_file', 'user_id', 'out_time', 'link')
    list_display_links = ('anal_name','insert_file', 'out_time')
    search_fields = ('anal_name','insert_file','user_id')
    list_per_page = 25
admin.site.register(ColapserDatabase, ColapserDatabaseAdmin)

class StrechingDatabaseAdmin(admin.ModelAdmin):
    list_display=('anal_name','user_id',)
    list_display_links = ('anal_name', 'user_id')
    search_fields = ('anal_name', 'user_id')
    list_per_page = 25 
admin.site.register(StretcherDatabase, StrechingDatabaseAdmin)
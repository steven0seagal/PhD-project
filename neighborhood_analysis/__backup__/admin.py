from django.contrib import admin

# Register your models here.
from .models import NeighAnalyzDatabase

class NeighAnalyzDatabaseAdmin(admin.ModelAdmin):
    list_display = ('id', 'user_id', 'pfam_name', 'neigh_size')
    list_display_links = ('id','user_id','pfam_name' )
    search_fields = ('user_id','pfam_name')
    list_per_page = 25
admin.site.register(NeighAnalyzDatabase, NeighAnalyzDatabaseAdmin)
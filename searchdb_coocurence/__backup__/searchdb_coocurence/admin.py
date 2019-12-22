from django.contrib import admin

# Register your models here.
from .models import Coocurence

class CoocurenceAdmin(admin.ModelAdmin):
    list_display = ('id', 'gene1', 'gene2')
    list_display_links = ('id','gene1', 'gene2' ) 
    
    list_per_page = 100
    search_fields = ('gene1', 'gene2', 'pvalue')
admin.site.register(Coocurence, CoocurenceAdmin )
from django.contrib import admin
from .models import HmmerFixerDatabase
# Register your models here.

class HmmerFixerDatabaseAdmin(admin.ModelAdmin):
    list_display = ('anal_name','main_file', 'user_id', 'skipped_seqs', 'link')
    list_display_links = ('anal_name','main_file', 'skipped_seqs')
    search_fields = ('anal_name','main_file','user_id','skipped_seqs')
    list_per_page = 25
admin.site.register(HmmerFixerDatabase, HmmerFixerDatabaseAdmin)
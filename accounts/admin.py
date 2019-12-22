from django.contrib import admin

# Register your models here.
from .models import PeekUserData, CompleteQueue

class PeekUserDataAdmin(admin.ModelAdmin):
    list_display = ('username', 'email')
    list_display_links = ('username', 'email')
    search_fields = ('username', 'email')
    list_per_page = 25


class CompleteQueueAdmin(admin.ModelAdmin):
    list_display = ('user_id','tool','status','file')
    list_display_links = ('user_id','tool','file')
    list_per_page = 25
    search_fields = ('user_id','status','tool','file')


admin.site.register(PeekUserData, PeekUserDataAdmin)
admin.site.register(CompleteQueue, CompleteQueueAdmin)
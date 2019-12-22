from django.db import models

# Create your models here.
class HmmerFixerDatabase(models.Model):
    main_file = models.CharField(blank=True, max_length = 200)
    skipped_seqs = models.CharField(blank=True, max_length = 200)
    user_id = models.IntegerField()
    anal_name = models.CharField(blank=True, max_length = 200)
    link = models.CharField(max_length = 200)
    def __str__(self):
        return self.anal_name
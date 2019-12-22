from django.db import models

# Create your models here.
class ColapserDatabase(models.Model):
    insert_file = models.CharField(max_length = 200)
    user_id = models.IntegerField()
    out_time = models.CharField(max_length = 50, blank=True)
    link = models.CharField(max_length = 200)
    anal_name = models.CharField(max_length = 200,blank = True)
    def __str__(self):
        return self.insert_file

class StretcherDatabase(models.Model):
    insert_file_master_master = models.CharField(max_length = 200)
    insert_file_master_slave_one = models.CharField(max_length = 200)
    insert_file_master_slave_two = models.CharField(max_length = 200)
    user_id = models.IntegerField()
    out_time = models.CharField(max_length = 50, blank=True)
    link = models.CharField(max_length = 200)
    anal_name = models.CharField(max_length = 200,blank = True)
    def __str__(self):
        return self.insert_file_master_master
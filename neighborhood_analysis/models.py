from django.db import models

# Create your models here.
class NeighAnalyzDatabase(models.Model):
    pfam_name = models.CharField(max_length = 9)
    neigh_size = models.IntegerField()
    cut_off = models.CharField(max_length = 9)
    tax = models.CharField(max_length = 200)
    strand = models.CharField(max_length = 200)
    out_name = models.CharField(blank=True, max_length = 200)
    user_id = models.IntegerField()
    out_time = models.CharField(max_length = 50, blank=True)
    link = models.CharField(max_length = 200)
    test_correction = models.CharField(max_length = 200)
    def __str__(self):
        return self.out_name
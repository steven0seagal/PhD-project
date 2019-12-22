from django.db import models

# Create your models here.

class Coocurence(models.Model):
    gene1 = models.CharField(max_length = 20 )
    gene2 = models.CharField(max_length = 20 )
    together  = models.IntegerField()
    first_only  = models.IntegerField()
    second_only  = models.IntegerField()
    neither  = models.IntegerField()
    pvalue  = models.DecimalField(max_digits=500, decimal_places=100)

def __str__(self):
    return self.gene1, self.gene2

from django.db import models

# Create your models here.

class Coocurence(models.Model):
    gene1 = models.CharField(max_length = 50 )
    gene2 = models.CharField(max_length = 50 )
    together  = models.IntegerField()
    first_only  = models.IntegerField()
    second_only  = models.IntegerField()
    neither  = models.IntegerField()
    pvalue  = models.DecimalField(max_digits=500, decimal_places=100)
    data_base = models.CharField(max_length = 100, blank = True )

class ColapsedOnSpeciesLevel(models.Model):
    gene1 = models.CharField(max_length = 50 )
    gene2 = models.CharField(max_length = 50 )
    together  = models.IntegerField()
    first_only  = models.IntegerField()
    second_only  = models.IntegerField()
    neither  = models.IntegerField()
    pvalue  = models.DecimalField(max_digits=500, decimal_places=100)
    data_base = models.CharField(max_length = 100, blank = True )

class ColapsedOnLegionellaStrains(models.Model):
    gene1 = models.CharField(max_length = 50 )
    gene2 = models.CharField(max_length = 50 )
    together  = models.IntegerField()
    first_only  = models.IntegerField()
    second_only  = models.IntegerField()
    neither  = models.IntegerField()
    pvalue  = models.DecimalField(max_digits=500, decimal_places=100)
    data_base = models.CharField(max_length = 100, blank = True )

class ColapsedOnLegionellaStrainWithingSpecies(models.Model):
    gene1 = models.CharField(max_length = 50 )
    gene2 = models.CharField(max_length = 50 )
    together  = models.IntegerField()
    first_only  = models.IntegerField()
    second_only  = models.IntegerField()
    neither  = models.IntegerField()
    pvalue  = models.DecimalField(max_digits=500, decimal_places=100)
    data_base = models.CharField(max_length = 100, blank = True )


def __str__(self):
    return self.gene1, self.gene2


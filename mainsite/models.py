
from django.db import models

class Tools(models.Model):
    image = models.ImageField(upload_to='photos/%Y/%m/%d/')
    summary = models.CharField(max_length = 200)
    strona = models.CharField(max_length = 50, default='SOME URL')

    def __str__(self):
        return self.summary
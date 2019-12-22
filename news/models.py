from django.db import models
from datetime import datetime
# Create your models here.
class NewsRoom(models.Model):
    author = models.CharField(max_length = 200 )
    title = models.CharField(max_length=100)
    pub_date = models.DateTimeField(default = datetime.now, blank = True )
    body = models.TextField()
    photo = models.ImageField(upload_to='photos/%Y/%m/%d/')
    is_published = models.BooleanField(default = True)
    para1 = models.TextField(blank = True)
    photo1 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para2 = models.TextField(blank = True)
    photo2 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para3 = models.TextField(blank = True)
    photo3 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para4 = models.TextField(blank = True)
    photo4 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para5 = models.TextField(blank = True)
    photo5 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para6= models.TextField(blank = True)
    photo6 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para7 = models.TextField(blank = True)
    photo7 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para8 = models.TextField(blank = True)
    photo8 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para9 = models.TextField(blank = True)
    photo9 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)
    para10 = models.TextField(blank = True)
    photo10 = models.ImageField(upload_to='photos/%Y/%m/%d/', blank = True)



    def summary(self):
        return self.body[:100]
    def pub_time_better(self):
        return self.pub_date.strftime('%a %d %B %Y')
    def __str__(self):
        return self.title
    def full_pub_time(self):
        return self.pub_date.strftime('%A %B %d %X')
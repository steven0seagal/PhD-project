from django.db import models


# Create your models here.
class PeekUserData(models.Model):


    first_name = models.CharField(max_length = 100)
    last_name = models.CharField(max_length = 100)
    username = models.CharField(max_length = 100)
    email = models.CharField(max_length = 100)
    password = models.CharField(max_length = 100)

    def __str__(self):
        return self.username


class CompleteQueue(models.Model):

    user_id = models.IntegerField()
    tool = models.CharField(max_length = 10)
    status = models.CharField(max_length = 25)
    analysis_name = models.CharField(max_length = 100)
    script = models.CharField(max_length = 500)
    file = models.CharField(max_length=200, blank= True)
    def __str__(self):
        return self.script

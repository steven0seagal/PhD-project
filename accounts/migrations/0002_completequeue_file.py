# Generated by Django 2.1.7 on 2019-09-24 12:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('accounts', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='completequeue',
            name='file',
            field=models.CharField(blank=True, max_length=200),
        ),
    ]
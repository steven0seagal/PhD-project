# Generated by Django 2.1.7 on 2019-04-19 10:29

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('searchdb_coocurence', '0005_auto_20190315_1900'),
    ]

    operations = [
        migrations.AddField(
            model_name='coocurence',
            name='data_base',
            field=models.CharField(blank=True, max_length=100),
        ),
    ]
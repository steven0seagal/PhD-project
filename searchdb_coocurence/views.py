from django.shortcuts import render
from django.http import HttpResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from searchdb_coocurence.choices import cutoff_choices, colapsing_choices
from .models import Coocurence,ColapsedOnSpeciesLevel,ColapsedOnLegionellaStrains,ColapsedOnLegionellaStrainWithingSpecies

def searchdbcooc(request):

    context = {
        'cutoff_choices':cutoff_choices,
        'colapsing_choices': colapsing_choices
    }

    return render(request, 'tools/input/coocurence.html',context)

def searchresult(request):

    colapsing  = request.GET['colapsing']

    if colapsing == '':
        queryset_list = Coocurence.objects.order_by('pvalue')
        
        if 'gene1' in request.GET:
            gene1 = request.GET['gene1']
            if gene1:
                queryset_list = queryset_list.filter(gene1__icontains = gene1)
        elif 'gene1' not in request.GET:
            gene1 = ''

        if 'gene2' in request.GET:
            gene2 = request.GET['gene2']
            if gene2:
                queryset_list = queryset_list.filter(gene2__icontains = gene2)
        else:
            gene2=''

        if 'pvalue' in request.GET:
            pvalue = request.GET['pvalue']
            if pvalue:
                queryset_list = queryset_list.filter(pvalue__lte=pvalue)
        else:
            pvalue= ''


    elif colapsing == 'Collapsed on species level':

        queryset_list = ColapsedOnSpeciesLevel.objects.order_by('pvalue')
        
        if 'gene1' in request.GET:
            gene1 = request.GET['gene1']
            if gene1:
                queryset_list = queryset_list.filter(gene1__icontains = gene1)
        elif 'gene1' not in request.GET:
            gene1 = ''

        if 'gene2' in request.GET:
            gene2 = request.GET['gene2']
            if gene2:
                queryset_list = queryset_list.filter(gene2__icontains = gene2)
        else:
            gene2=''

        if 'pvalue' in request.GET:
            pvalue = request.GET['pvalue']
            if pvalue:
                queryset_list = queryset_list.filter(pvalue__lte=pvalue)
        else:
            pvalue= ''

    elif colapsing == 'Collapsed on Legionella strains':

        queryset_list = ColapsedOnLegionellaStrains.objects.order_by('pvalue')
        
        if 'gene1' in request.GET:
            gene1 = request.GET['gene1']
            if gene1:
                queryset_list = queryset_list.filter(gene1__icontains = gene1)
        elif 'gene1' not in request.GET:
            gene1 = ''

        if 'gene2' in request.GET:
            gene2 = request.GET['gene2']
            if gene2:
                queryset_list = queryset_list.filter(gene2__icontains = gene2)
        else:
            gene2=''

        if 'pvalue' in request.GET:
            pvalue = request.GET['pvalue']
            if pvalue:
                queryset_list = queryset_list.filter(pvalue__lte=pvalue)
        else:
            pvalue= ''

    elif colapsing == 'Collapsed on species within Legionella':

        queryset_list = ColapsedOnLegionellaStrainWithingSpecies.objects.order_by('pvalue')
        
        if 'gene1' in request.GET:
            gene1 = request.GET['gene1']
            if gene1:
                queryset_list = queryset_list.filter(gene1__icontains = gene1)
        elif 'gene1' not in request.GET:
            gene1 = ''

        if 'gene2' in request.GET:
            gene2 = request.GET['gene2']
            if gene2:
                queryset_list = queryset_list.filter(gene2__icontains = gene2)
        else:
            gene2=''

        if 'pvalue' in request.GET:
            pvalue = request.GET['pvalue']
            if pvalue:
                queryset_list = queryset_list.filter(pvalue__lte=pvalue)
        else:
            pvalue= ''





    page = request.GET.get('page',1)
    paginator = Paginator(queryset_list, 100)
    
    query = paginator.get_page(page)

    context = {
         'gene1':gene1,
         'gene2':gene2,
         'pvalue':pvalue,
         'colapsing':colapsing,
         'query':query,
         'cutoff_choices':cutoff_choices,
         'colapsing_choices': colapsing_choices,
         
     }
    return render(request, 'tools/output/coocurence_result.html' ,context)


    #return render(request, 'tools/output/coocurence_result.html', context)

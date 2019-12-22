from django.shortcuts import render
from django.http import HttpResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from searchdb_coocurence.choices import cutoff_choices, taxa_choices
from .models import Coocurence
def searchdbcooc(request):

    context = {
        'cutoff_choices':cutoff_choices,
    }

    return render(request, 'tools/input/coocurence.html',context)

def searchresult(request):

    queryset_list = Coocurence.objects.order_by('pvalue')

    if 'gene1' in request.GET:
        gene1 = request.GET['gene1']
        if gene1:
            queryset_list = queryset_list.filter(gene1__icontains = gene1)

    if 'gene2' in request.GET:
        gene2 = request.GET['gene2']
        if gene2:
            queryset_list = queryset_list.filter(gene2__icontains = gene2)

    if 'pvalue' in request.GET:
            pvalue = request.GET['pvalue']
            if pvalue:
                queryset_list = queryset_list.filter(pvalue__lte=pvalue)




    paginator = Paginator(queryset_list, 5000)
    page = request.GET.get('page')
    query = paginator.get_page(page)

    context = {
         'query':query,
         'cutoff_choices':cutoff_choices,
     }
    return render(request, 'tools/output/coocurence_result.html' ,context)


    #return render(request, 'tools/output/coocurence_result.html', context)

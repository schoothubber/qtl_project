from django.conf.urls import patterns,url
from qtl import views


urlpatterns = patterns('',
	url(r'^trait/$', views.SearchTraitView, name = 'search_individual_trait'),
	url(r'^check/$', views.MultipleTraitView, name = 'search_multiple_trait'),
	url(r'^graph/$', views.Graphview, name='display_graph'),
	url(r'^godag/$', views.Graphview, name='display_godag'),
	url(r'^output/$', views.OutputDataView, name='display_data'),
	url(r'^storecsv/$', views.StoreCsvView, name='store_csv'),
						) 

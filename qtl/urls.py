from django.conf.urls import patterns,url
from qtl import views


urlpatterns = patterns('',
	url(r'^trait/$', views.SearchTraitView, name = 'search_trait'),
	url(r'^check/$', views.CheckTraitView, name = 'check_trait'),
						) 

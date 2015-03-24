from django.conf.urls import patterns, include, url
from qtl.views import BaseView

urlpatterns = patterns('',
	url(r'^$', BaseView, name='home'),
	url(r'^eqtlagogo/', include('qtl.urls')),
	)
	

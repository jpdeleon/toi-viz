from django.urls import path
#added for media
from django.conf import settings
from django.conf.urls.static import static

from . import views

#add namespace to differentiate apps in the project
app_name = 'polls'
urlpatterns = [
    # ex: /polls/
    #path('', views.index, name='index'),
    path('', views.IndexView.as_view(), name='index'),
    # ex: /polls/specifics/5/
    #path('specifics/<int:question_id>/', views.detail, name='detail'),
    path('<int:pk>/', views.DetailView.as_view(), name='detail'),
    # ex: /polls/5/results/
    #path('<int:question_id>/results/', views.results, name='results'),
    path('<int:pk>/results/', views.ResultsView.as_view(), name='results'),
    # ex: /polls/5/vote/
    #path('<int:question_id>/vote/', views.vote, name='vote'),
    path('<int:question_id>/vote/', views.vote, name='vote'),
] + static(settings.MEDIA_URL, document_root= settings.MEDIA_ROOT) # added for media


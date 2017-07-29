from django.views.generic import FormView

from . import forms


# Create your views here.
class UploadView(FormView):
    template_name = 'upload/index.html'
    form_class = forms.UploadForm
    success_url = 'http://coruzzilab-macpro.bio.nyu.edu'

    def form_valid(self, form: forms.UploadForm):
        form.save_files()
        form.send_mail()
        return super().form_valid(form)

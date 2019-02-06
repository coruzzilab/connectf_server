import gzip
import json
import os
import shutil
from glob import iglob

import pandas as pd
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.core.files.storage import FileSystemStorage
from django.test import TestCase
from django.urls import reverse

from querytgdb.utils.insert_data import import_additional_edges, import_annotations, insert_data, read_annotation_file
from .models import Analysis, Annotation, EdgeData, EdgeType


# Create your tests here.
class TestImportData(TestCase):
    @classmethod
    def setUpClass(cls):
        annotation_file = max(iglob("test_data/annotation*.csv.gz"), key=os.path.getmtime)
        import_annotations(annotation_file)

    @classmethod
    def tearDownClass(cls):
        Annotation.objects.all().delete()

    def test_import_annotations(self):
        Annotation.objects.all().delete()

        annotation_file = max(iglob("test_data/annotation*.csv.gz"), key=os.path.getmtime)

        annotations = read_annotation_file(annotation_file)

        import_annotations(annotation_file)

        new_annotations = pd.DataFrame(Annotation.objects.values_list(named=True).iterator()).drop('id', axis=1)

        self.assertTrue(annotations.equals(new_annotations), "newly imported data should be the same as file")

    def test_import_additional_edges(self):
        with gzip.open('test_data/dap_tga1.txt.gz', 'r') as g:
            import_additional_edges(g, sif=False, directional=True)

        with self.subTest("should have new edges"):
            self.assertTrue(EdgeType.objects.filter(name='DAP').exists(), "should have DAP edge type")
            self.assertTrue(EdgeType.objects.filter(name='ampDAP').exists(), "shoud have ampDAP edge type")
            self.assertTrue(EdgeData.objects.filter(tf__gene_id="AT5G65210").exists(),
                            "should have new TGA1 additional edges")

    def test_import_experiment_data(self):
        with self.assertRaises(ObjectDoesNotExist):
            Analysis.objects.get(tf__gene_id="AT5G65210")

        with gzip.open('test_data/AT5G65210_TGA1_DESEQ2.txt.gz', 'rt') as m, \
                gzip.open('test_data/AT5G65210_TGA1_DESeq2.csv.gz', 'r') as d:
            insert_data(d, m)

        with self.subTest("should have new analysis"):
            analysis = Analysis.objects.get(tf__gene_id="AT5G65210")

            self.assertTrue(analysis.analysisdata_set.exists(), "should have metadata")
            self.assertTrue(analysis.interaction_set.exists(), "should have target data")
            self.assertTrue(analysis.regulation_set.exists(), "should have p-value and fold change data")


static_storage = FileSystemStorage(settings.QUERY_CACHE)


class TestQuery(TestCase):
    @classmethod
    def setUpClass(cls):
        annotation_file = max(iglob("test_data/annotation*.csv.gz"), key=os.path.getmtime)
        import_annotations(annotation_file)

        with gzip.open('test_data/AT5G65210_TGA1_DESEQ2.txt.gz', 'rt') as m, \
                gzip.open('test_data/AT5G65210_TGA1_DESeq2.csv.gz', 'r') as d:
            insert_data(d, m)

        with gzip.open('test_data/dap_tga1.txt.gz', 'r') as g:
            import_additional_edges(g, sif=False, directional=True)

    @classmethod
    def tearDownClass(cls):
        pass

    def test_basic_query(self):
        response = self.client.post(reverse("queryapp:queryapp"), data={
            "query": "AT5G65210"
        })

        result = json.loads(response.content)

        self.assertIn('request_id', result, "has request_id")

        cache_path = static_storage.path(f"{result['request_id']}_pickle")

        self.assertTrue(os.path.exists(cache_path))

        with self.subTest("all files have data"):
            cache_files = [
                'formatted_tabular_output.pickle.gz',
                'metadata.pickle.gz',
                'query.txt',
                'tabular_output.pickle.gz',
                'tabular_output_unfiltered.pickle.gz',
            ]
            for f in cache_files:
                self.assertGreater(os.path.getsize(os.path.join(cache_path, f)), 0, f"{f} should have data")

        shutil.rmtree(cache_path)

    def test_cached_query(self):
        response = self.client.post(reverse("queryapp:queryapp"), data={
            "query": "AT5G65210"
        })

        request_id = json.loads(response.content)['request_id']

        response = self.client.get(reverse("queryapp:queryapp") + request_id + "/")

        self.assertEqual(response.status_code, 200)

        shutil.rmtree(static_storage.path(f"{request_id}_pickle"))

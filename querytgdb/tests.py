import gzip
import io
import json
import os
import secrets
from glob import iglob

import pandas as pd
from django.core.exceptions import ObjectDoesNotExist
from django.test import TestCase
from django.urls import reverse

from querytgdb.utils.insert_data import import_additional_edges, import_annotations, insert_data, \
    read_annotation_file
from .models import Analysis, Annotation, EdgeData, EdgeType
from .utils.file import BadNetwork, get_network


class TestImportData(TestCase):
    @classmethod
    def setUpClass(cls):
        annotation_file = os.environ.get('TEST_ANNOTATION',
                                         max(iglob("test_data/annotation*.csv.gz"), key=os.path.getmtime))
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
        """
        Check if all cache files exist with basic query
        :return:
        """
        response = self.client.post(reverse("queryapp:queryapp"), data={
            "query": "AT5G65210"
        })

        result = json.loads(response.content)

        self.assertIn('request_id', result, "has request_id")

    def test_cached_query(self):
        """
        Test if cached query is used
        :return:
        """
        response = self.client.post(reverse("queryapp:queryapp"), data={
            "query": "AT5G65210"
        })

        request_id = json.loads(response.content)['request_id']

        response = self.client.get(reverse("queryapp:queryapp") + request_id + "/")

        self.assertEqual(response.status_code, 200)


class TestNetworkParsing(TestCase):
    def test_good_file(self):
        buff = io.StringIO("source	DFG_Prediction	dest	score\n"
                           "AT4G25210	DFG_Prediction	AT4G13940	54.252\n"
                           "AT4G36540	DFG_Prediction	AT4G13940	44.818")
        name, data = get_network(buff)

        self.assertEqual(name, "default", "Should have default name")
        self.assertIsInstance(data, pd.DataFrame, "Should be dataframe")
        self.assertEqual(data.shape, (2, 5), "should have 2 rows 5 columns")

    def test_bad_file(self):
        buff = io.BytesIO(secrets.token_bytes(1024))  # if this turns out to be a valid network, go buy a lottery ticket

        with self.assertRaises(BadNetwork):
            get_network(buff)

    def test_empty_file(self):
        buff = io.BytesIO()

        with self.assertRaises(BadNetwork):
            get_network(buff)

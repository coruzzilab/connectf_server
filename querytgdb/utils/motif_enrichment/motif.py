import re
from abc import ABC
from collections import OrderedDict
from typing import Dict, List, Optional, Type

import pandas as pd
import seaborn as sns
from django.conf import settings

from querytgdb.utils import async_loader, skip_for_management


class MotifError(Exception):
    pass


@skip_for_management
def get_annotations():
    return pd.read_csv(settings.MOTIF_ANNOTATION, index_col=[0, 1, 2], header=None)


@skip_for_management
def get_tf_annotations():
    return pd.read_csv(settings.MOTIF_TF_ANNOTATION, index_col=[0, 1, 2], header=None)


async_loader['motifs'] = get_annotations
async_loader['motifs_tf'] = get_tf_annotations


class MotifData:
    _regions: Dict[str, 'Region'] = OrderedDict()

    def __init__(self, background: Optional[pd.Series] = None):
        self.cache: Dict[str, pd.DataFrame] = {}
        self._colors = None
        self.background = background

        self._annotation = None

    def get_annotation(self) -> pd.DataFrame:
        """
        Override for alternative source of motif annotations
        :return:
        """
        return async_loader['motifs']

    @property
    def annotation(self) -> pd.DataFrame:
        if self._annotation is not None:
            return self._annotation

        anno = self.get_annotation()

        if self.background is not None:  # account for background
            anno = anno[anno.index.get_level_values(0).isin(self.background)]

        self._annotation = anno

        return self._annotation

    def __getattr__(self, item):
        try:
            return self.cache[item]
        except KeyError:
            if item.endswith('_cluster_size'):
                return self.cluster_size(re.sub(r'_cluster_size$', '', item))

            return self.get_region(item)

    def __getitem__(self, item) -> 'Region':
        return self._regions[item]

    def cluster_size(self, region):
        try:
            return self.cache[region + '_cluster_size']
        except KeyError:
            cluster_size = getattr(self, region).groupby(level=2).sum()
            self.cache[region + '_cluster_size'] = cluster_size

            return cluster_size

    @property
    def region_total(self):
        try:
            return self.cache['region_total']
        except KeyError:
            region_total = self.annotation.groupby(level=1).sum()
            self.cache['region_total'] = region_total
            return region_total

    def get_region(self, name):
        region_matches = self.annotation.loc[(slice(None), self._regions[name].name, slice(None)), :]
        self.cache[name] = region_matches

        return region_matches

    @classmethod
    def register(cls, region_cls: Type['Region']) -> Type['Region']:
        region = region_cls()

        name = region.name

        cls._regions[name] = region

        return region_cls

    @property
    def motifs(self) -> pd.Index:
        return self.annotation.index.levels[2]

    @property
    def regions(self) -> List[str]:
        return list(self._regions.keys())

    @property
    def colors(self):
        if self._colors is None:
            self._colors = dict(zip(self.regions, sns.color_palette("husl", len(self.regions))))

        return self._colors

    @property
    def default_regions(self) -> List[str]:
        return [key for key, val in self._regions.items() if val.default]

    @property
    def default_region(self) -> str:
        try:
            return self.default_regions[0]
        except IndexError:
            raise MotifError("No default region")

    @property
    def region_desc(self) -> Dict[str, Dict]:
        """
        Include discription of each region
        :return:
        """
        return OrderedDict((key, val.to_dict()) for key, val in self._regions.items())


class AdditionalMotifData(MotifData):
    def get_annotation(self) -> pd.DataFrame:
        return async_loader['motifs_tf']


class Region(ABC):
    default: bool = False
    name: str = ''
    description: str = ''
    group: List = []

    def __init__(self):
        if not self.name:
            self.name = self.__class__.__name__.lower()

    def __repr__(self):
        return f'<Region: {self.name}>'

    def to_dict(self) -> Dict:
        return {
            'name': self.name,
            'description': self.description,
            'group': self.group
        }

import pathlib
import re
from abc import ABC, abstractmethod
from collections import OrderedDict
from threading import Lock, Thread
from typing import Dict, List, Type, Union

import pandas as pd


class MotifData:
    def __init__(self, data_file: Union[str, pathlib.Path]):
        self.lock = Lock()
        self.data_file = data_file
        self._annotation = None

        self.task = Thread(target=self.get_annotation)
        self.task.start()

        self._regions: Dict[str, 'Region'] = OrderedDict()

        self.cache: Dict[str, pd.DataFrame] = {}

    def get_annotation(self):
        self._annotation = pd.read_csv(self.data_file, index_col=0)

    @property
    def annotation(self) -> pd.DataFrame:
        with self.lock:
            if self._annotation is None:
                self.task.join()

        return self._annotation

    def __getattr__(self, item):
        try:
            return self.cache[item]
        except KeyError:
            if item.endswith('_dedup'):
                return self.dedup(re.sub(r'_dedup$', '', item))
            elif item.endswith('_cluster_size'):
                return self.cluster_size(re.sub(r'_cluster_size$', '', item))

        raise AttributeError(item)

    def __getitem__(self, item) -> 'Region':
        return self._regions[item]

    def dedup(self, region):
        try:
            return self.cache[region + '_dedup']
        except KeyError:
            dedup = getattr(self, region).drop_duplicates('match_id')
            self.cache[region + '_dedup'] = dedup
            return dedup

    def cluster_size(self, region):
        try:
            return self.cache[region + '_cluster_size']
        except KeyError:
            cluster_size = self.dedup(region).groupby('#pattern name').size()
            self.cache[region + '_cluster_size'] = cluster_size

            return cluster_size

    def register(self, region_cls: Type['Region']) -> Type['Region']:
        region = region_cls()

        name = region.name

        self._regions[name] = region
        self.cache[name] = region.get_region(self.annotation)

        return region_cls

    @property
    def regions(self) -> List[str]:
        return list(self._regions.keys())

    @property
    def default_regions(self) -> List[str]:
        return [key for key, val in self._regions.items() if val.default]

    @property
    def region_desc(self) -> Dict[str, Dict]:
        """
        Include discription of each region
        :return:
        """
        return OrderedDict((key, val.to_dict()) for key, val in self._regions.items())


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

    @abstractmethod
    def get_region(self, annotation: pd.DataFrame) -> pd.DataFrame:
        raise NotImplementedError("Should implement a region filter.")

    def to_dict(self) -> Dict:
        return {
            'name': self.name,
            'description': self.description,
            'group': self.group
        }

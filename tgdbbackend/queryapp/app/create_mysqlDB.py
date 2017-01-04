#!/usr/bin/python
# -*- coding: utf-8 -*-

''' This script setup the mysql database for Target system DB
'''

############
# Modules
############

import sys
import argparse as agp
from sqlalchemy import Column, ForeignKey, Integer, String, alias, Text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
from warnings import filterwarnings

Base = declarative_base()
filterwarnings('ignore')

class Nodes(Base):

	__tablename__ = 'Nodes'
	node_id = Column(Integer, primary_key=True, autoincrement=True, index= True)
	node_name = Column(String(100), unique=True, nullable=False, index= True)
	node_type = Column(String(100), nullable=False, index= True)

class Edges(Base):

	__tablename__= 'Edges'
	edge_id= Column(Integer, primary_key=True, autoincrement=True, index= True)
	edge_name= Column(String(100), unique=True, nullable=False, index= True)

class Meta(Base):

	__tablename__='Meta'
	meta_p= Column(Integer, primary_key= True, autoincrement=True) # I had to make this a primary key
	meta_id= Column(String(100), primary_key= True, nullable=False, index= True)
	analysis_id= Column(String(100), primary_key= True, nullable=False, index= True)
	meta_value= Column(String(255),index=True) # meta value should be a text
	meta_type= Column(String(100),nullable=False, index= True)

class Interactions(Base):

	__tablename__ = 'Interactions'
	interaction_id = Column(Integer, primary_key=True, autoincrement=True)
	node_1_id = Column(Integer, ForeignKey('Nodes.node_id'), nullable=False, index= True)
	node_2_id = Column(Integer, ForeignKey('Nodes.node_id'), nullable=False, index= True)
	edge_id = Column(Integer, ForeignKey('Edges.edge_id'), nullable=False)
	meta_id= Column(String(100), ForeignKey('Meta.meta_id'), nullable=False) # could not change this to an interger because of autoincrement
	analysis_id= Column(String(100), ForeignKey('Meta.analysis_id'), nullable=False)

class Genenames(Base):
	
	__tablename__ = 'Genenames'
	ath_id= Column(String(100), primary_key=True, unique=True, index= True)
	ath_name= Column(String(200))
	ath_fullname= Column(Text(2000))
	ath_gene_type= Column(String(100),nullable=False)
	ath_gene_fam= Column(Text(2000),nullable=False)

class createdbase():

	def __init__(self,dbasename):
		self.engine = create_engine('mysql+pymysql://root:coruzzilab@localhost')
		self.engine.execute("DROP DATABASE IF EXISTS "+dbasename) # uncomment this if you want to recreate the database
		self.engine.execute("CREATE DATABASE "+dbasename)
		self.engine.execute("USE "+dbasename)
		Base.metadata.create_all(self.engine)
		
	def getEngine(self):
		return self.engine

class AccessDatabase(): # This class has special permissions for coruzzilab users(have only permissions to access the database)

	def __init__(self, dbname):
		self.engine = create_engine('mysql+pymysql://coruzzilab:accesstargetdb@172.22.2.137/'+dbname)
	
	def getEngine(self):
		return self.engine

if __name__=='__main__':

	parser= agp.ArgumentParser()
	parser.add_argument('-d',help='Assign name to the mysql database')
	args = parser.parse_args()
	createdbase(args.d)


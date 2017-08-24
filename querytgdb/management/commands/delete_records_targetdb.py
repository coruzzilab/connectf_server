#!/usr/bin/python

'''
This script is to discard records (experiments or analysis) from the targetdb database.
-e option: If only experiment id is provided, it will delete the whole experiment (means all the analysis in an exp).
-a option: If analysis id is provided along with the experiment id, it will delete the individual analysis.


'''

import os,sys
import argparse as agp
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, aliased
from create_mysqlDB import Nodes,Edges,Meta,Interactions,Genenames,Base,createdbase,AccessDatabase


# delete an analysis
def del_analysis(sess, experimentid, analysisid):

    sess.query(Interactions).filter(Interactions.meta_id == experimentid). \
        filter(Interactions.analysis_id== analysisid). \
        delete(synchronize_session='fetch')
    sess.query(Meta).filter(Meta.meta_id == experimentid). \
        delete(synchronize_session='fetch')


# delete an experiment
def del_exp(sess, experimentid):
	print 'experimentid= ',experimentid

	#sess.query(Interactions).filter(Interactions.meta_id == experimentid).\
	#							delete(synchronize_session='fetch')

	sess.query(Meta).filter(Meta.meta_id == experimentid).\
								delete(synchronize_session='fetch')


############################
# main method
#@profile
def main(dbname, experiment_id, analysis_id):

	engine = create_engine('mysql://root:coruzzilab@localhost/'+dbname)
	Base.metadata.bind= engine
	DBSession= sessionmaker(bind=engine)
	sess= DBSession()


	if experiment_id and analysis_id: # delete an analysis
		del_analysis(sess, experiment_id, analysis_id)



	if not analysis_id: # delete the whole experiment
		del_exp(sess, experiment_id)


	sess.commit()



if __name__=='__main__':

	parser= agp.ArgumentParser()
	parser.add_argument('-d','--dbname',help= 'Database name',required= True)
	# program only deletes one experiment/analysis at a time, so the argument provided are not lists
	parser.add_argument('-e','--experimentid', help= 'delete an experiment',required= True)
	parser.add_argument('-a','--analysisid', help= 'delete an analysis within an experiment')

	args= parser.parse_args()

	main(args.dbname, args.experimentid, args.analysisid)



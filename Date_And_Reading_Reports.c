
/********************************************
*	Date_And_Reading_Reports.c				*
*	Created by Neal Smith in March 2018.	*
*	Used to convert a date DD/MM/YYYY into	*
*		a date_ID, starting from some		*
*		specific date.						*
********************************************/

//preprocessor directives
#include <stdio.h>						//For standard input/output functions
#include <stdlib.h>						//For standard C library functions
#include <string.h>						//For functions that manipulate strings
#include "Date_And_Reading_Reports.h"	//For structures and declarations of functions needed in this file
											//To allow multiple source files to be used in one program
#include "MTrandom.h"					//For random number generation (accept/reject situations)
#include "lfunc.h"						//For MTrandom.cpp

/*---------------------------------------------------------------
| functions contained in this source code, for use in Ebola_x.c |
---------------------------------------------------------------*/

int generate_report_date_id(int day, int month, int year)
{
	int value = 0;

	//This date ID is counted from the likely infectious period start date of the index case (Dec 2, 2013?)
	value += 365 * (year - 2014);
	value += 30 * month;
	if (month >= 1) {
		if (month >= 2) {
			if (month >= 3) {
				if (month >= 4) {
					if (month >= 6) {
						if (month >= 8) {
							if (month >= 9) {
								if (month >= 11) {
									value += 1;
								}
								value += 1;
							}
							value += 1;
						}
						value += 1;
					}
					value += 1;
				}
				value -= 2;
			}
			value += 1;
		}
		value += 1;
	}
	value += day - 2 + 1;

	return value;
}

int assign_x() //F:: NEEDS PROPER INPUT VARIABLES
{
	printf("Assigning x coordinate for case based on report location.\n");
		//This will be removed when proper function is written	
	return 0;	//Temporary - until we can assign a proper value
}

int assign_y()	//F:: NEEDS PROPER INPUT VARIABLES
{
	printf("Assigning y coordinate for case based on report location.\n");
	//This will be removed when proper function is written	
	return 0;	//Temporary - until we can assign a proper value
}

void generate_cases(struct current_case_report *current_report, p_patient *head, struct parameter_list *p_params)
{	//I've passed the above to this function as it is outside the main source code and these parts are needed
	int i;						//to count the number of cases created so far
	int j;						//to count days when needed to initialise values
	p_patient current;	//the current case we're working with
	p_patient index_case;	//the first case of the outbreak

	if (p_params->total_cases == 0){	//if it's the first case of all,
		//creating index case
		*head = index_case = (p_patient*)malloc(sizeof(p_patient));
		if(!index_case){printf("Could not allocate index case in generate_cases.\n"); exit(1);}

		//Initialising index case variables
		strcpy(index_case->country, current_report->country);		//Country of case
		strcpy(index_case->subregion, current_report->subregion);	//Subregion of case
		//index_case->x = assign_x();
		//index_case->x = (int)(uniform() * 100);		//Gives case an x-coord b/w 1 + 100 (0 and 99?)
		//index_case->y = assign_y();
		//index_case->y = (int)(uniform() * 100);		//Gives case an x-coord b/w 1 + 100 (0 and 99?)
		index_case->diag = 1;		//This case is from a report, so diagnosis = 1
		index_case->diag_day = 1;	//Index case is diagnosed on date_ID = 1
		index_case->est_case = 1;	//Any reported case became established
		index_case->parent_case = NULL;	//No parent case for index case
		index_case->secondary_cases = 0;	//No secondary cases yet
		for (j = 1; j <= NUMDAYS; j++) {
			index_case->secondary_cases_gen[j] = 0;	//No secondary cases on any day
		}
		index_case->first_2dary = NULL;	//No cases founded
		index_case->left_sib = NULL;		//No cases founded
		index_case->right_sib = NULL;		//No cases founded
		index_case->prev = NULL;			//No-one newer in the list
		index_case->next = NULL;			//No-one older in the list
		index_case->transmission_type = 0;//Zoonotic transmission for index case (0)
		//index_case->pop_dens = population_density(x, y)	//Likely format for population density input
			//VARIABLES NOT DEALT WITH: DATES[N], SURVIVAL
	}
	else {	//creating non-index case
		current = (struct patient*)malloc(sizeof(struct patient));
		if (!current){printf("Could not allocate new case in generate_cases.\n"); exit(1);}
		
		//Initialising index case variables
		strcpy(current->country, current_report->country);		//Country of case
		strcpy(current->subregion, current_report->subregion);	//Subregion of case
		//current->x = assign_x();
		//current->x = (int)(uniform()*100);		//Gives case an x-coord b/w 1 + 100 (0 and 99?)
		//current->y = assign_y();
		//current->y = (int)(uniform()*100);		//Gives case an x-coord b/w 1 + 100 (0 and 99?)
		current->diag = 1;		//This case is from a report, so diagnosis = 1
		//current->diag_day = assign_diag_day();	//From date_ID and prev_date_ID (for country?)
		current->est_case = 1;	//Any reported case became established
		current->secondary_cases = 0;	//No secondary cases yet
		for (j = 1; j <= NUMDAYS; j++) {
			current->secondary_cases_gen[j] = 0;	//No secondary cases on any day
		}
		current->first_2dary = NULL;	//No cases founded
		current->left_sib = NULL;		//No cases founded
		current->right_sib = NULL;		//No cases founded
		current->transmission_type = 0;//Zoonotic transmission for index case (0)
		//current->pop_dens = population_density(x, y)	//Likely format for population density input
			//VARIABLES NOT DEALT WITH: DATES[N], SURVIVAL, PARENT
			
			//Generate linked list
		current->next = *head;		//next in the list is the most recently created case
		current->prev = NULL;		//no previous case as this is the most recent
		(*head)->prev = current;		//if head is next, current is head's previous
		*head = current;				//switch head to most current (to become new case's next case) 
	}
		(p_params->num_diagnosed)++;//New case - update diagnosed case count
		(p_params->total_cases)++;	//New case - update total case count
}
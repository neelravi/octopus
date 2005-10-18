/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id$
*/

#include <config.h>

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "string_f.h" /* fortran <-> c string compatibility issues */

typedef struct opt_type{
  char *name;
  char *value;
  char *desc;
  struct opt_type *next;
} opt_type;

typedef struct var_type{
  char *name;
  char *type;
  char *section;
  char *desc;
  opt_type *opt;
  struct var_type *next;
} var_type;

static var_type *vars = NULL;


/* --------------------------------------------------------- */
char *get_token(char *s, char **dest)
{
  char *s1;

  /* get rid of initial whitespace */
  for(;*s!='\0' && isspace(*s); s++); 
  if(!isalnum(*s) && *s!='-'){
    *dest = NULL;
    return s;
  }

  for(s1=s; isalnum(*s1) || *s1=='_' || *s1=='-'; s1++);

  *dest = strndup(s, s1-s);
  return s1;
}


/* --------------------------------------------------------- */
void get_text(FILE *in, char **dest)
{
  char c, line[256];

  for(;;){
    /* check if the next line starts by a space */
    if((c=getc(in)) == EOF) return;
    ungetc(c, in);

    if(!isspace(c)) return;

    fgets(line, 256, in);
    if(!*dest)
      *dest = strdup(line+1);
    else{
      *dest = realloc(*dest, strlen(*dest)+strlen(line+1)+1); 
      strcat(*dest, line+1);
    }
  }
}


/* --------------------------------------------------------- */
void FC_FUNC_(varinfo_init, VARINFO_INIT)
  (STR_F_TYPE fname STR_ARG1)
{
  char line[256], *fname_c;
  FILE *in;
  var_type *lvar = NULL;
  opt_type *lopt;

  fname_c = TO_C_STR1(fname);

  in = fopen(fname_c, "r");
  if(!in) goto out;

  while(fgets(line, 256, in)){

    if(strncasecmp("Variable", line, 8) == 0){
      char *s;

      get_token(line+9, &s);
      if(s){ /* found a token */
	if(!lvar){
	  lvar = (var_type *) malloc(sizeof(var_type));
	  vars = lvar;
	}else{
	  lvar->next = (var_type *) malloc(sizeof(var_type));
	  lvar = lvar->next;
	}
	lvar->name = s;
	lvar->desc = NULL;
	lvar->type = NULL;
	lvar->section = NULL;
	lvar->opt  = NULL;
	lvar->next = NULL;

	lopt = NULL;
      }
      continue;
    }

    /* if no variable was found continue */
    if(!lvar) continue;

    if(strncasecmp("Type", line, 4) == 0)
      get_token(line+5, &(lvar->type));

    if(strncasecmp("Section", line, 7) == 0){
      char *s = line+7;
      for(; *s!='\0' && isspace(*s); s++);
      lvar->section = strdup(s);
    }

    if(strncasecmp("Description", line, 11) == 0){
      if(lvar->desc){ /* if repeated delete old description */
	free(lvar->desc);
	lvar->desc = NULL;
      }
      get_text(in, &(lvar->desc));
    }

    if(strncasecmp("Option", line, 6) == 0){
      char *name, *value, *s;
      s = get_token(line+6, &name);
      if(name) get_token(s, &value);

      if(name){ /* found an option */
	if(!lopt){
	  lopt = (opt_type *) malloc(sizeof(opt_type));
	  lvar->opt = lopt;
	}else{
	  lopt->next = (opt_type *) malloc(sizeof(var_type));
	  lopt = lopt->next;
	}
	lopt->name  = name;
	lopt->value = value;
	lopt->desc  = NULL;
	get_text(in, &(lopt->desc));
	lopt->next  = NULL;	
      }
    }
  }
  fclose(in);
  
 out:
  free(fname_c);
}


/* --------------------------------------------------------- */
void FC_FUNC_(varinfo_end, VARINFO_END)
  ()
{
  var_type *v = vars;
  for(;v;){
    var_type *v1 = v->next;
    opt_type *o  = v->opt;

    if(v->name) free(v->name);
    if(v->type) free(v->type);
    if(v->section) free(v->section);
    if(v->desc) free(v->desc);
    for(;o;){
      opt_type *o1 = o->next;
      if(o->name ) free(o->name);
      if(o->value) free(o->value);
      if(o->desc)  free(o->desc);

      free(o);
      o = o1;
    }

    free(v);
    v = v1;
  }
}


/* --------------------------------------------------------- */
void FC_FUNC_(varinfo_getvar, VARINFO_GETVAR)
  (STR_F_TYPE name, var_type **var STR_ARG1)
{
  char *name_c;
  var_type *lvar;

  name_c = TO_C_STR1(name);
  for(lvar=vars; (lvar!=NULL) && (strcasecmp(name_c, lvar->name)!=0); lvar=lvar->next);
  free(name_c);

  *var = lvar;
}


/* --------------------------------------------------------- */
void FC_FUNC_(varinfo_getinfo, VARINFO_GETVAR)
  (var_type **var, char **name, char **type, char **section, char **desc)
{
  if(var == NULL){
    *name = NULL; *type = NULL; *desc = NULL;
  }else{
    *name = (*var)->name;
    *type = (*var)->type;
    *section = (*var)->section;
    *desc = (*var)->desc;
  }
}


/* --------------------------------------------------------- */
void FC_FUNC_(varinfo_getopt, VARINFO_GETOPT)
  (var_type **var, opt_type **opt)
{
  if(*var == NULL)
    *opt = NULL;
  else if(*opt == NULL)
    *opt = (*var)->opt;
  else
    *opt = (*opt)->next;
}


/* --------------------------------------------------------- */
void FC_FUNC_(varinfo_opt_getinfo, VARINFO_OPT_GETINFO)
  (opt_type **opt, char **name, int *value, char **desc)
{
  if(opt == NULL){
    *name = NULL; *desc = NULL;
    *value = 0;
  }else{
    *name = (*opt)->name;
    *desc = (*opt)->desc;
    if((*opt)->value)
      *value = atoi((*opt)->value);
    else
      *value = 0;
  }
}

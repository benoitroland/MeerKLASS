# lotss-query

This repository contains the module that enables interaction with the
LoTSS surveys database. ddf-pipeline depends on this module but it can
be used stand alone.

Users must have an account on the LoTSS database server (or run their
own local equivalent). By default, without local environment variables
set, the code will try to establish an ssh tunnel to the LoTSS
database server and use that for database communication.

# Database documentation

The LoTSS HBA database back end is a MySQL database with 6 tables:

```
MariaDB [surveys]> show tables;
+-------------------+
| Tables_in_surveys |
+-------------------+
| fields            |
| observations      |
| quality           |
| reprocessing      |
| spectra           |
| transients        |
+-------------------+
```

These respectively contain:

* fields: the properties and current status of each of the LoTSS fields
* observations: observations taken. Each observation must be of one of the fields.
* quality: the output of the quality pipeline.
* reprocessing: lists the current reprocessing tasks and their status.
* spectra: the current list of dynamical spectra and their properties (each dynamical spectrum should correspond to a source in the transients table)
* transients: the current list of sources for which dynamical spectra should be found.

You can access other databases by the `survey` optional argument to `SurveysDB`, e.g. `survey='lba'`. 

Programmers' notes
------------------

Make sure that your home directory contains a file called `.surveys`
with the database password, and optionally your remote username and
ssh key in it.

All the work to provide the python interface is in the surveys_db module.

```
from surveys_db import SurveysDB

with SurveysDB() as sdb:
    result=sdb.db_get('fields','P35Hetdex10')
    print(result)
```

The SurveysDB class does the work of setting up an ssh tunnel (if run
remotely) and creating a connection to the database. Note that it
locks the two tables for writing (by default: change with
'readonly=True') so only one instance can exist anywhere at any
time. So you should remember to use a context manager (as above) or
run the `close()` method as soon as you're finished with it.  If using
the ssh tunnel method a fixed port is used by default so only one
instance per machine is possible (change with localport=xxxxx).

SurveysDB exports the `.cur` method which is a PyMySQLdb cursor, so
low-level operations are possible. If you want to carry out a complex query, you should design it in MySQL and use the `.cur` method.

Or you can use the `db_get`, `db_set` and `db_create` methods. These take two arguments,
a table to use and a dictionary or an ID. `db_get` populates a
dictionary with the results of a query for that table and ID: dictionary keys
are column names. `db_set` takes a dictionary and puts it back into
the corresponding table: the dictionary key `id` must be an existing ID and all
other keys must be valid columns. And `db_create` makes a new table
row and returns a dictionary which can be updated and passed to
`db_set`.

Methods `get_field`, `set_field` etc are wrappers around
`db_get('field', ...)` etc and should not be used in new code.

Use in the pipeline
-------------------

To use the database hooks in ddf-pipeline itself set the environment variable `DDF_PIPELINE_DATABASE`.

Setting this variable has the following effects:

* download_field script will register a download of a field.
* unpack script will register unpacking
* make_mslist script will register that the data are ready
* pipeline will update status at the start and end of the run and on failure

You should also set `DDF_PIPELINE_CLUSTER` to the name of your cluster, e.g. 'Herts', 'Leiden'.


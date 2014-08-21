drop table if exists projects;

create table projects (
    id integer primary key autoincrement,
    title text not null,
    model text not null,
    N1 integer not null,
    max_snps integer not null,
    max_mstats integer not null,
    replicates integer not null,
    max_cycles integer not null,
    burn_in integer not null)

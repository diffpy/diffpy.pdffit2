// Numbering of the parameters

// n_st    : Number of global structural parameters per phase
// n_at    : Number of parameters for each atom
// n_ex    : Number of experimental parameters per dataset

phase.offset
atom[ia].offset = phase.offset + n_st + ia*n_at;
set[is].offset

       parameter       (n_ex =  3)
        parameter       (n_st = 10)
        parameter       (n_at = 10)

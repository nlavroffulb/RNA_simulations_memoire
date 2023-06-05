#pragma once
#include <vector>
//double p_select_helix_origin() { return 0.25; }

double p_select_linked_structure();

double p_select_unlinked_structure();

double Z_ev_linked();
double Z_ev_unlinked();

double partition_ni_unlinked(int length);

double partition_ni_linked(int length);

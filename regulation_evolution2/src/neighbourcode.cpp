/*
This code should be somewhere with access to the CA plane and the cells --either dish or ca
*/

/*Call this function after placing cells in the plane, if you want neighbour info*/
void Agent::InitContactLength() {
    int i, j, k;
    int celltype, celltypeneigh;
    int boundary = 0;
    int totalperimeter = 0;
    int sigma, sigmaneigh;
    map<int, Cell>::iterator it;

    for (it = CellArray.begin(); it != CellArray.end(); ++it) {
        it->second.clearNeighbours();
    }

    for (i = 0; i < L; i++) {
        for (j = 0; j < W; j++) {
            if ((sigma = CellIdGrid[i][j]) && (celltype = CellArray[CellIdGrid[i][j]].celltype)) //focus is on a cell
            {
                for (k = 0; k < neighbourhoodsize; k++)//go through neighbourhood of the pixel
                {
                    if ((sigmaneigh = CellIdGrid[i + wideneighbourhood[k].yy][j + wideneighbourhood[k].xx]) !=
                        sigma)//medium can also be a neighbour!
                    {
                        totalperimeter++;
                        CellArray[sigma].updateNeighbour(sigmaneigh, 1);
                    }
                }
            }
        }
    }

}

//after divisions, neighbours need to be re-established. !!! need to have the old position of the cell
// so call this function before UpdateCell in the division functions.
void Agent::RecountContactLength(int id, int newid) {
    int i, j, k;
    int sigma, sigmaneigh;

    if (id == 0 || newid == 0) {
        printf("RecountContactLength: error, updating medium for contacts\n");
        exit(1);
    }

    set<int> store;
    set<int> store2;
    double scale = 3.;
    ///find the section of the field over which to iterate (taken quite big)
    int minL = (int) max(CellArray[id].meani - scale * CellArray[id].ReturnMajorAxisLength(), 4.);
    int maxL = (int) min(CellArray[id].meani + scale * CellArray[id].ReturnMajorAxisLength(), (double) L - 5);
    int minW = (int) max(CellArray[id].meanj - scale * CellArray[id].ReturnMajorAxisLength(), 4.);
    int maxW = (int) min(CellArray[id].meanj + scale * CellArray[id].ReturnMajorAxisLength(), (double) W - 5);

    ///go through the neighbours of the old cell, and remove the contact with ID
    map<int, int>::iterator it;
    int checkbounds = 0;
    int count;
    for (it = CellArray[id].neighbours.begin(), count = 0; it != CellArray[id].neighbours.end(); ++it, count++) {
        checkbounds += it->second;
        store.insert(it->first);
        if (it->first)
            CellArray[it->first].setNeighbour(id, 0); //this neighbour removes cell ID from its contacts

    }
    CellArray[id].clearNeighbours();

    int check = 0;
    int checkbounds2 = 0;
    ///reset the contact by iterating over the field
    for (i = minL; i <= maxL; i++) {
        for (j = minW; j <= maxW; j++) {
            if ((sigma = CellIdGrid[i][j]) == id || sigma == newid) //focus is on the cell or its daughter
            {
                for (k = 0; k < neighbourhoodsize; k++)//go through neighbourhood of the pixel
                {
                    if ((sigmaneigh = CellIdGrid[i + wideneighbourhood[k].yy][j + wideneighbourhood[k].xx]) != sigma) {
                        check += CellArray[sigma].updateNeighbour(sigmaneigh, 1);
                        if (sigmaneigh != id && sigmaneigh != newid) {
                            checkbounds2++;
                            store2.insert(sigmaneigh);
                        }
                        if (sigmaneigh != id && sigmaneigh != newid &&
                            sigmaneigh) {//also update the contacts of the neighbour, unless it is one of the two focus cells:otherwise you update those twice
                            check += CellArray[sigmaneigh].updateNeighbour(sigma, 1);
                        }
                    }
                    if (check) {
                        printf("error in RecountContactlength: updating pixels sigma %d and sigmaneigh %d\n", sigma,
                               sigmaneigh);
                        printf("agent nr %d\n", agentid);
                    }
                }
                if (i == minL || i == maxL || j == minW || j == maxW)
                    printf("RecountContactLength: caution: cell in question at border of field. cell: %d, pos: %d %d\n",
                           id, i, j);
            }
        }
    }

    if (checkbounds2 != checkbounds) {
        printf("error in RecountContactLength: wrong nr of boundary pixels! divided %d to %d,  checkbounds %d, checkbounds2 %d\n",
               id, newid, checkbounds, checkbounds2);
        for (i = minL; i < maxL; i++) {
            for (j = minW; j < maxW; j++) {
                printf("%d\t", CellIdGrid[i][j]);
            }
            printf("\n");
        }
        // exit(1);
    }
    if (store != store2) {
        printf("error in RecountContactLength: different registered neighbours!\n");
        ///for debugging:
//     set<int> diff1;
//     printf("in store 1 but not store 2: ");
//     set_difference(store.begin(), store.end(), store2.begin(), store2.end(),inserter(diff1, diff1.begin()));
//   
//     for (set<int>::iterator dit=diff1.begin(); dit!=diff1.end(); ++dit) printf("%d\t", (*dit));
//     printf("\n");
//     
//     set<int> diff2;
//     printf("in store 2 but not store 1: ");
//     set_difference(store2.begin(), store2.end(), store.begin(), store.end(),inserter(diff2, diff2.begin()));
//   
//     for (set<int>::iterator dit=diff2.begin(); dit!=diff2.end(); ++dit) printf("%d\t", (*dit));
//     printf("\n");
        InitContactLength();
    }

}


/**************************************/
/*
This code should be in the function that does the bookkeeping after a pixel copy
*/
/**************************************/

int check = 0;
///use this if you need info about contacts with neighbours
for(
int k = 0;
k<neighbourhoodsize;
k++)
{
point = CellIdGrid[i + wideneighbourhood[k].yy][j + wideneighbourhood[k].xx];

if(point!=idA)
{
if(
idA &&CellArray
.
count(idA)
)
check+=CellArray[idA].
updateNeighbour(point,
-1);
if(point)
check+=CellArray[point].
updateNeighbour(idA,
-1);
}
if(point!=idB)
{
if(
idB &&CellArray
.
count(idB)
)
check+=CellArray[idB].
updateNeighbour(point,
1);
if(point)
check+=CellArray[point].
updateNeighbour(idB,
1);
}
if(check)
{
printf("error in After: wrongly updating neighbours of copy event idA %d and idB %d\n", idA, idB);
printf("agent nr %d\n", agentid);
}
}



/**************************************/
/*
This code should be in the cell header file
*/
/**************************************/

map<int, int> neighbours; //stores neighbouring cells(ID) and the amount of membrane contact
//cell neighbours
void setNeighbour(int neighbour, int boundarylength);

void clearNeighbours();

int returnBoundaryLength(int cell);

int updateNeighbour(int cell, int modification);



/**************************************/
/*
This code should be in the cell cpp file
*/
/**************************************/


void Cell::setNeighbour(int neighbour, int boundarylength) {

    if (boundarylength == 0)//remove this neighbour
        neighbours.erase(neighbour);
    else
        neighbours[neighbour] = boundarylength; //if the element is already present, the boundarylength will be modified, otherwise a new element will be created.

}

int Cell::returnBoundaryLength(int cell) {
    if (neighbours.count(cell))
        return neighbours[cell];

    return 0;

}

void Cell::clearNeighbours() {
    neighbours.clear();
}

int Cell::updateNeighbour(int cell, int modification) {
    if (!neighbours.count(cell) && modification < 0) {
        printf("Cell.updateNeighbour: error: negatively updating contact of cell %d with nonexisting neighbour %d\n",
               id, cell);
        return 1;
    } else if (!neighbours.count(cell))
        neighbours[cell] = modification;

    else if (neighbours.count(cell)) {
        neighbours[cell] += modification;
        if (neighbours[cell] == 0)//remove this neighbour
        {
            neighbours.erase(cell);
        } else if (neighbours[cell] < 0) {
            neighbours.erase(cell);
            printf("Cell.updateNeighbour: error: updating contact of cell %d with neighbour %d to negative value\n", id,
                   cell);
            return 2;
        }
    }
    return 0;
}





#include "formangradient.h"

void FormanGradient::computeBoundaryCell(implicitS const& cell, SSet& descCells){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    descCells = SSet(foo);

    stack<implicitS> qu;
    qu.push(cell);

    while(!qu.empty()){

        implicitS top = qu.top();
        qu.pop();

        vector<implicitS>* bd = sc.boundaryk(top,top.getDim()-1);
        for(auto s : *bd){
            implicitS next;
            if(getPair(s,next)){
                assert(isPaired(next) && isPaired(s));
                if(next.getDim() == cell.getDim() && !sc.theSame(next,top)){
                    qu.push(next);
                }
            }
            else{
                if(!isPaired(s)){
                    if(descCells.find(s) == descCells.end())
                        descCells.insert(s);
                    else
                        descCells.erase(s);
                }
            }
        }
        delete bd;
    }
}

#include "helpers.hpp"
#include "structures.hpp"

//helper functions

void group_bodies(Rect rect, Body * bodies, int n, 
                  Rect * rects, Body ** groups, int * ns){
        //init rects
        Real minx = rect.pos0.x;
        Real miny = rect.pos0.y;
        Real maxx = rect.pos1.x;
        Real maxy = rect.pos1.y;
        Real midx = (minx+maxx)/2;
        Real midy = (miny+maxy)/2;
        
        rects[0] = {rect.pos0,{midx,midy}};
        rects[1] = {{midx,miny},{maxx,midy}};
        rects[2] = {{minx,midy},{midx,maxy}};
        rects[3] = {{midx,midy},rect.pos1};
        //recursively call
        for(int i = 0; i < n; i++){
                for(int j = 0; j < 4; j++){
                        if(rects[j].contains(bodies[i].p)){
                                groups[j][ns[j]] = bodies[i];
                                ns[j] ++;
                                break;
                        }
                }
        }
}

Tree built_qtree_r(Rect rect, Body * bodies, int n){
        Tree tree_children[4];
        Node node = {.rect = rect};
        node.body = NULL;
        Tree tree = {.node = node, .trees = tree_children};
        Rect rects[4];
        Body groups[4][n];
        int ns[4];
        group_bodies(rect, bodies, n, rects, groups, ns);
        //recursively call
        int num_child = 4;
        for(int j = 0; j < 4; j++){
                if(ns[j] > 1){
                        tree.trees[j] = build_qtree_r(rects[j], groups[j], ns[j]);
                } else if (ns[j]){
                        Body body = groups[j][0];
                        tree.trees[j] = {{rects[j], {body.p, 0, body.mass}, body},NULL};
                } else {
                        num_child --;
                }
        }
        //compute center
        tree.num_child = num_child;
        Rvec * centers = Rvec[num_child];
        for(int j = 0; j < num_child; j++){
                centers[j] = tree.trees[j].center;
        }
        tree.center = mix_rvecs(centers, num_child);
        //return
        return tree;
}

void NBody::metric_expansion(Real dt){
        Real ratio = 1+hubble*dt;
        uargs.size *= ratio;
        for(int i = 0; i < num_bodies; i++){
                bodies[i].p = bodies[i].p.M(ratio);
        }
}

void NBody::border_wrap(void){
        for(int i = 0; i < num_bodies; i++){
                bodies[i].p.x = fmod(bodies[i].p.x, uargs.size);
                bodies[i].p.y = fmod(bodies[i].p.y, uargs.size);
        }
}

void NBody::build_qtree(void){
        return build_qtree_r({{0,0},{size,size}}, bodies, uargs->numbodies);
}

void NBody::leapfrog(Real dt){
        Vec * accs = accelerations();
        for(int i = 0; i < num_bodies; i++){
                bodies[i].p.P(bodies[i].v.M(dt));
                bodies[i].v.P(accs[i].M(dt));
        }
}

Vec NBody::accel(Vec p1, Vec p2, Real mass, Real distance){
        Real mag_accel = uargs.gconst*mass/(distance*distance+plummer*plummer);
        Real dir_accel = atan2(p2.y-p1.y, p2.x-p1.x);
        return {acc*cos(dir_accel), acc*sin(dir_accel)};
}

Vec NBody::lattice_accel(Vec pos, RVec center, Real lattice_dist){
        Vec accel = {0,0};
        int lat = uargs->lattice;
        //find all csom and distances
        Vec * CsOM = Vec[lat*lat];
        Real * distances = Real[lat*lat];
        for(int x = -lat; i <= lat; i++){
                for(y = -lat; y <= lat; y++){
                        int index = x*lat+y;
                        Vec Com = {pos.x+x*size, pos.y+y*size};
                        CsOM[index] = COM;
                        distances[index] = distance(pos, COM);
                }
        }
        //compute acceleration
        for(int i = 0; i < lat*lat; i++){
                //why?
                if(distances[i]/lattice_dist < lattice+1){
                        accel.P(pos, center.vec, center.mass, distances[i]);
                }
        }
        return accel;
}

Vec NBody::body_accel(int id){
        Vec accel = {0,0};
        static Tree_stack BFSS = {0, Tree[uargs->nbody]};
        BFSS.add(quadtree);
        while(BFSS.n){
                tree = BFSS.pop();
                Real lattice_dist = lat_dist(bodies[i].p, tree.center.var, uargs->size);
                if(tree.node.body != NULL && tree.node.body->id == id){
                        continue;
                } else if (tree.node.center.var == 0 || 
                           lattice_dist/sqrt(tree.center.var) > QTR){
                        accel.P(lattice_accel(bodies[i].p, tree.node.center, lattice_dist));
                } else {
                        for(int i = 0; i < tree.num_child; i++){
                                BFSS.add(tree.trees[i]);
                        }
                }
        }
        return accel;
}

Vec * NBody::accelerations(void){ 
        static Vec * accs;
        if(accs == NULL){
                Vec * accs = new Vec[num_bodies];
        }
        for(int i = 0; i < num_bodies; i++){
                accs[i] = body_accel(i);
        }
        return accs;
}

void NBody::step(Real dt, Dargs dargs = NULL){
        build_qtree();
        leapfrog(dt);
        border_wrap();
        metric_expansion(dt);
        if(dargs != NULL){
                draw(dargs);
        }
}

void NBody::simulate(Real simtime, Real dt, string simname, int disp_size, verbose = True){
        for(int i = 0; i < simetime/dt; i++){
                if(verbose){
                        printf("%f", dt*i);
                }
                step(dt, {simname+to_string(i), disp_size});
        }
}

void NBody::draw(Dargs dargs){
        Real dsize = dargs.size;
        static bool * parr;
        if(parr == NULL){
                parr = new bool[dsize][dsize];
        }
        memset(parr, false, dsize*dsize*sizeof(bool));
        
        for(int i = 0; i < uargs.num_bodies; i++){
                int x = bodies[i].p.x/uargs.size*dsize;
                int y = bodies[i].p.x/uargs.size*dsize;
                to_image[x%dsize][y%dsize] = true;
        }
        to_image(parr, dsize, dargs.filename);
}
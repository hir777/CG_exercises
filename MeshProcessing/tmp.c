int i, num_tris;
    getNumberTriangles(&g_tri_mesh, &num_tris);

    glShadeModel(GL_SMOOTH);
    glBegin(GL_TRIANGLES);
    for( i = 0; i < num_tris; i++ ) {
      Vector3 v[3], n;
      getTriangleVertices(&g_tri_mesh, i, v);

      getVertexNormal(&g_tri_mesh, g_tri_mesh._triangles[i]._v0, &n);
      glNormal3f(n._x, n._y, n._z);
      glVertex3f(v[0]._x, v[0]._y, v[0]._z);

      getVertexNormal(&g_tri_mesh, g_tri_mesh._triangles[i]._v1, &n);
      glNormal3f(n._x, n._y, n._z);
      glVertex3f(v[1]._x, v[1]._y, v[1]._z);

      getVertexNormal(&g_tri_mesh, g_tri_mesh._triangles[i]._v2, &n);
      glNormal3f(n._x, n._y, n._z);
      glVertex3f(v[2]._x, v[2]._y, v[2]._z);
    }
    glEnd();




int i, j;
  float weight = 1.0;

  int num_verts = tri_mesh->_number_vertices;
  int num_tris  = tri_mesh->_number_triangles;
  tri_mesh->_vertex_normals = (Vector3 *)malloc( num_verts * sizeof(Vector3) )
  Vector3 sum;
  for( i = 0; i < num_verts; i++ ) {
    sum._x = sum._y = sum._z = 0.0;
    Vector3 v = tri_mesh->_vertices[i];

    for( j = 0; j < num_tris; j++ ) {
      Vector3 coord[3];
      getTriangleVertices(tri_mesh, j, coord);
      if( compare(v, coord[0]) || compare(v, coord[1]) || compare(v, coord[2]) ) {
        // true only when the triangle is adjacent to the vertex v
        Vector3 n;
        getTriangleNormal(tri_mesh, j, &n);
        mulAV(weight, n, &n);
        add(sum, n, &sum);
      }
    }

    normalize(sum, &sum);
    tri_mesh->_vertex_normals[i]._x = sum._x;
    tri_mesh->_vertex_normals[i]._y = sum._y;
    tri_mesh->_vertex_normals[i]._z = sum._z;
  }




// The function compares the two given vector.
// Return value:
// 1 equal
// 0 not-equal
int compare(Vector3 v1, Vector3 v2) {
  if( v1._x == v2._x && v1._y == v2._y && v1._z == v2._z )
    return 1;
  else
    return 0;
}

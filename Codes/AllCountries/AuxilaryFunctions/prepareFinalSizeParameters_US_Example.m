   f   u   n   c   t   i   o   n       [   u   n   i   f   o   r   m   A   l   l   o   c   a   t   i   o   n   ,   r   i   s   k   F   o   c   u   s   e   d   A   l   l   o   c   a   t   i   o   n   ,   s   p   r   e   a   d   e   r   s   A   l   l   o   c   a   t   i   o   n   ,   u   p   p   e   r   B   o   u   n   d   ,   v   a   c   c   i   n   e   s   L   e   f   t   T   o   D   i   s   t   i   b   u   t   e   ,   a   d   u   l   t   A   g   e   s   ,   a   g   e   A   b   o   v   e   6   0   ,   a   g   e   4   0   t   o   6   0   ,   a   g   e   2   0   t   o   4   0   ,   I   F   R   ,   C   i   j   ,   N   i   ,   N   a   d   u   l   t   ,   r   ,   v   ,   i   n   f   e   c   t   e   d   0   _   v   ,   i   n   f   e   c   t   e   d   0   _   n   v   ]   =   p   r   e   p   a   r   e   F   i   n   a   l   S   i   z   e   P   a   r   a   m   e   t   e   r   s   _   U   S   _   E   x   a   m   p   l   e   (   V   c   P   r   c   t   C   h   i   l   d   r   e   n   ,   b   e   t   a   V   a   c   ,   e   f   f   V   a   c   ,   m   a   x   P   r   c   t   )   
   %       c   o   u   n   t   r   y   L   i   s   t   =   {   "   B   E   L   "   ,       "   U   S   A   "   ,       "   I   N   D   "   ,       "   E   S   P   "   ,       "   Z   W   E   "   ,       "   B   R   A   "   ,       "   C   H   N   "   ,       "   Z   A   F   "   ,       "   P   O   L   "   }   ;   
   %       R   0   =   3   ;   
   %       r   e   c   o   v   e   r   e   d   p   r   c   t   =   0   ;   
   %       i   n   f   e   c   t   e   d   _   n   v   _   p   r   c   t   =   0   ;   
   %       i   n   f   e   c   t   e   d   _   v   _   p   r   c   t   =   0   ;   
   %       V   c   P   r   c   t   =   8   0   ;       %       T   h   e       p   e   r   c   e   n   t       o   f       2   0   +       p   o   p   u   l   a   t   i   o   n       t   h   a   t       h   a   v   e       v   a   c   c   i   n   e   s       a   v   a   l   i   a   b   l   e   
   %       b   e   t   a   V   a   c   =   0   .   2   ;       %       R   e   l   a   t   i   v   e       s   u   s   c   e   b   t   i   b   i   l   i   t   y       o   f       v   a   c   c   i   n   e   d       i   n   d   i   v   i   d   u   a   l   s   
   %       e   f   f   V   a   c   =   0   .   9   5   ;       %       E   f   f   i   c   a   c   y       o   f       v   a   c   c   i   n   e       i   n       p   r   e   v   e   n   t   i   n   g       d   i   s   e   a   s   e   
   %       m   a   x   P   r   c   t   =   9   5   ;       %       M   a   x   i   m   a   l       v   a   c   c   i   n   a   t   i   o   n       p   e   r       a   g   e       g   r   o   u   p   
   
   %   %       L   o   a   d       c   o   u   n   t   r   y       d   a   t   a   
   c   o   u   n   t   r   y   D   a   t   a   =   l   o   a   d   (   j   o   i   n   (   [   '   .   .   /   c   o   u   n   t   r   y   D   a   t   a   /   U   S   A   _   d   a   t   a   .   m   a   t   '   ]   ,   '   '   )   )   ;   
   
   %   %       P   r   e   p   a   r   e       d   a   t   a   
   C   i   j   =   c   o   u   n   t   r   y   D   a   t   a   .   c   o   n   t   a   c   t   M   a   t   r   i   x   ;   
   I   F   R   =   c   o   u   n   t   r   y   D   a   t   a   .   I   F   R   '   ;   
   N   =   c   o   u   n   t   r   y   D   a   t   a   .   N   ;   
   N   i   =   N   *   c   o   u   n   t   r   y   D   a   t   a   .   a   g   D   i   s   t   ;   
   a   g   e   A   b   o   v   e   2   0   =   3   :   9   ;   a   g   e   A   b   o   v   e   6   0   =   7   :   9   ;   a   g   e   4   0   t   o   6   0   =   5   :   6   ;   a   g   e   2   0   t   o   4   0   =   3   :   4   ;   
   a   g   e   A   b   o   v   e   6   0   =   a   g   e   A   b   o   v   e   6   0   -   2   ;   a   g   e   4   0   t   o   6   0   =   a   g   e   4   0   t   o   6   0   -   2   ;   a   g   e   2   0   t   o   4   0   =   a   g   e   2   0   t   o   4   0   -   2   ;       %       i   n   d   e   x       1       r   e   f   e   r   s       t   o       f   i   r   s   t       a   d   u   l   t       a   g   e   
   
   u   p   p   e   r   b   o   u   n   d   F   i   l   t   e   r   =   o   n   e   s   (   s   i   z   e   (   N   i   )   )   ;   
   N   a   d   u   l   t   =   N   i   ;   
   a   d   u   l   t   A   g   e   s   =   1   :   9   ;   
   
   %       C   o   m   p   u   t   e       n   e   x   t       g   e   n   e   r   a   t   i   o   n       m   a   t   r   i   x       &       s   c   a   l   e       c   o   n   t   a   c   t       m   a   t   r   i   x   
   M   i   j   =   d   i   a   g   (   N   i   )   *   C   i   j   *   d   i   a   g   (   1   .   /   N   i   )   ;   
   [   V   ,   d   ]   =   e   i   g   (   M   i   j   )   ;   
   C   i   j   =   C   i   j   /   d   (   1   )   ;   
   i   n   f   e   c   t   e   d   D   i   s   t   r   i   b   u   t   i   o   n   =   V   (   :   ,   1   )   /   s   u   m   (   V   (   :   ,   1   )   )   ;       %       E   x   t   r   a   c   t       a   g   e       d   i   s   t   r   i   b   u   t   i   o   n       o   f       i   n   f   e   c   t   v   i   e   s       f   r   o   m       l   a   r   g   e   s   t       e   .   v   .       o   f       n   g   m   
   
   %       I   n   i   t   i   a   l   i   z   e       s   u   s   c   e   p   t   i   b   l   e   ,       r   e   c   o   v   e   r   e   d   ,       v   a   c   c   i   n   a   t   e   d       a   n   d       a   c   t   i   v   e       i   n   f   e   c   t   i   v   e   s   
   
   v   =   [   0   ;   0   .   6   2   ;   0   .   7   ;   0   .   7   ;   0   .   7   8   ;   0   .   8   6   ;   0   .   9   9   ;   0   .   9   8   ;   0   .   9   6   ]   ;   
   r   =   [   0   .   0   5   7   0   5   0   6   2   8   ;   0   .   1   0   9   1   6   3   8   7   1   ;   0   .   1   4   7   8   5   9   3   7   4   ;   0   .   1   3   7   9   7   1   6   5   8   ;   0   .   1   3   3   6   8   9   1   3   ;   0   .   1   2   0   5   2   1   3   6   8   ;   0   .   0   9   3   2   1   0   0   6   5   ;   0   .   0   8   0   9   9   4   5   1   9   ;   0   .   0   9   7   8   9   6   7   8   3   ]   
   V   c   =   [   0   .   5   *   N   i   (   1   )   +   N   i   (   2   )   *   0   .   2   ]   *   V   c   P   r   c   t   C   h   i   l   d   r   e   n   /   1   0   0   ;   
   s   =   m   a   x   (   1   -   r   -   v   ,   0   )   ;   
   
   %       C   o   m   p   u   t   e       u   n   i   f   o   r   m       a   l   l   o   c   a   t   i   o   n       o   f       v   a   c   c   i   n   e   s   
   a   g   e   F   i   l   t   e   r   =   0   *   N   i   ;   a   g   e   F   i   l   t   e   r   (   a   d   u   l   t   A   g   e   s   )   =   t   r   u   e   ;   
   v   a   c   c   i   n   e   s   L   e   f   t   T   o   D   i   s   t   i   b   u   t   e   =   V   c   ;   
   u   p   p   e   r   B   o   u   n   d   =   m   a   x   (   s   (   a   d   u   l   t   A   g   e   s   )   -   (   1   -   m   a   x   P   r   c   t   /   1   0   0   )   ,   0   )   ;   u   p   p   e   r   B   o   u   n   d   =   m   i   n   (   u   p   p   e   r   b   o   u   n   d   F   i   l   t   e   r   (   a   d   u   l   t   A   g   e   s   )   ,   u   p   p   e   r   B   o   u   n   d   )   ;   
   u   n   i   f   o   r   m   A   l   l   o   c   a   t   i   o   n   =   c   o   m   p   u   t   e   U   n   i   f   o   r   m   A   l   l   o   c   a   t   i   o   n   (   s   ,   a   d   u   l   t   A   g   e   s   ,   N   a   d   u   l   t   ,   v   a   c   c   i   n   e   s   L   e   f   t   T   o   D   i   s   t   i   b   u   t   e   ,   u   p   p   e   r   B   o   u   n   d   )   ;   
   r   i   s   k   F   o   c   u   s   e   d   A   l   l   o   c   a   t   i   o   n   =   c   o   m   p   u   t   e   R   i   s   k   F   o   c   u   s   e   d   A   l   l   o   c   a   t   i   o   n   (   s   ,   a   d   u   l   t   A   g   e   s   ,   N   a   d   u   l   t   ,   v   a   c   c   i   n   e   s   L   e   f   t   T   o   D   i   s   t   i   b   u   t   e   ,   u   p   p   e   r   B   o   u   n   d   )   ;   
   s   p   r   e   a   d   e   r   s   A   l   l   o   c   a   t   i   o   n   =   c   o   m   p   u   t   e   S   p   r   e   a   d   e   r   s   A   l   l   o   c   a   t   i   o   n   (   s   ,   a   d   u   l   t   A   g   e   s   ,   N   a   d   u   l   t   ,   v   a   c   c   i   n   e   s   L   e   f   t   T   o   D   i   s   t   i   b   u   t   e   ,   u   p   p   e   r   B   o   u   n   d   )   ;   
   i   n   f   e   c   t   e   d   0   _   v   =   0   *   N   i   ;   i   n   f   e   c   t   e   d   0   _   n   v   =   0   *   N   i   ;   
   r   e   t   u   r   n   

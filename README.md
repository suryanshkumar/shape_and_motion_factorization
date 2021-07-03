Statement: The above code implements the Tomasi et.al [1] classical factorization method to
rigid structure from motion under orthographic camera model assumption.

An initial version of this code was written by me at the start of my Ph.D. at 
Australian National University in October 2015. Improved and refined recently for 
better understanding. (Removed the Newton optimization module).

# code structure
1. dataset folder contains the .mat file of the tracked points. (Dataset downloaded from the internet.)
2. src folder contains the implementation of the algorithm.
3. run tomasi_kanade_factorization_1992.m file to see the results.

# Reference
<table>
<tr>
<td>
[1] <strong>Shape and Motion from Image Streams under Orthography: a Factorization Method.</strong><br />
Carlo Tomasi, Takeo Kanade <br /> International Journal of Computer Vision, IJCV, 1992, pages 137-154, Springer.<br />
[<a href="https://www-users.cs.umn.edu/~hspark/CSci5980/tomasi.pdf" target="_blank">pdf</a>] <br />
</td>
</tr>
</table>
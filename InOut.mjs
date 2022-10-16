// Function to process in out parameters

export function InOutF( /* in */ v1, /* in */ v2, /* in out */ io1, /* in out */ io2 )
{
    console.log( "INOUTF: In Out parameters BEFORE assignment" );
    console.log( "    v1 = ", v1 );
    console.log( "    v2 = ", v2 );
    console.log( "    io1 = ", io1 );
    console.log( "    io2 = ", io2 );
    // Perform assignment
    io1.Double = v1;
    io2.Double = v2;
    console.log( "INOUTF: In Out parameters AFTER assignment" );
    console.log( "    io1 = ", io1 );
    console.log( "    io2 = ", io2 );
}

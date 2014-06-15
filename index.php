<?
////////////////////////////////////////////////////////////
//ENTRADAS
////////////////////////////////////////////////////////////
foreach(array_keys($_GET) as $field){
    $$field=$_GET[$field];
}
foreach(array_keys($_POST) as $field){
    $$field=$_POST[$field];
}
$dut=-5;
?>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
</head>
<center>
<img src="images/aristarco.jpg" height=200>
<H1>Ocultación de Marte por la Luna<br/>
Julio 5 de 2014</H1>
</center>

<p>En esta página podrás calcular las condiciones de la ocultación de
Marte por la Luna del próximo 5 de Julio de 2014 para tu ubicación
específica.  Los observadores que se encuentran ubicados sobre la
franja de totalidad (ver mapa abajo) pueden determinar aquí los
tiempos de los contactos teniendo en cuenta o no el efecto de adelanto
del fenómeno debido a la velocidad finita de la luz.</p>

<H3><a name="Tiempos">Tiempos de Contacto</a></H3>

<p>A continuación ingrese su posición geográfica exacta (tal y como es
provista por un GPS).</p>
<form action="?#Tiempos">
Longitud: <input type="text" name="lon" value="<?echo $lon?>"><br/>
<i style=font-size:12px>Formato: +GG.gggg, +GG:MM:SS.sss, +GG MM SS.sss.  Ejemplos: -75 35 23.99, -75:35:23.99, -75.43555</i>
<p></p>
Latitud: <input type="text" name="lat" value="<?echo $lat?>"><br/>
<i style=font-size:12px>Formato: +GG.gggg, +GG:MM:SS.sss, +GG MM SS.sss.  Ejemplos: -75 35 23.99, -75:35:23.99, -75.43555</i>
<p></p>
Uso horario: <input type="text" name="dut" value="<?echo $dut?>"><br/>
<i style=font-size:12px>Formato: +UTC.  Ejemplo: -5, -4.5</i>
<p></p>
Altura: <input type="text" name="alt" value="<?echo $alt?>"> metros<br/>
<i style=font-size:12px>Formato: +AAAA.aaa.  No use "," ni separadores de miles.</i>
<p></p>
<input type='submit' name='accion' value='Calcule'>
</form>
<?
////////////////////////////////////////
//CALCULO
////////////////////////////////////////
if($accion=="Calcule")
{
  //PREPARE INPUTS
  
  //CALL PROGRAM
  $out=shell_exec("./contact-times.out $lon $lat $alt $dut 2> tmp/contact-times.log");

  //GENERATE CORD
  $filename=sprintf("data/cord-%+.4lf_%+.4lf_%+.3lf.dat",$lon,$lat,$alt);
  shell_exec("python plot-cords.py $filename");
}
?>
<H3>Franja de Totalidad</H3>
<center><img src="images/totalidad.png"></center>
</html>

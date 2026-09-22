# Configurar rclone para este proyecto

Los scripts `drive_push.sh` y `drive_pull.sh` mueven cientos de GB entre la
máquina local y `Mi unidad/tesis/`. Esta guía es lo que ellos necesitan, no una
guía genérica de rclone.

**Cuando termines, corré `./scripts/drive_check.sh`.** Verifica la config contra
lo que los scripts esperan y distingue los dos modos de fallo que *no dan
error*. Si dice `TODO OK`, está.

---

## Lo que los scripts exigen, exactamente

| | valor | de dónde sale |
| :-- | :-- | :-- |
| nombre del remoto | **`gdrive-tesis`** | `DRIVE_REMOTE_DEFAULT` en `scripts/_drive_lib.sh` |
| tipo | **`drive`** | — |
| cuenta | **`seb.ugazm@gmail.com`** | `data/DRIVE.md` |
| carpeta raíz | **`tesis`** (como *ruta*, no como `root_folder_id`) | `DRIVE_ROOT` |

Con eso, un `.sra` de `prupe` vive en `gdrive-tesis:tesis/80_sra/prupe/`.

Los dos primeros se pueden cambiar sin tocar código: `DRIVE_REMOTE=otro` y
`DRIVE_ROOT=otra_carpeta` son variables de entorno.

---

## 1. Instalar rclone

```bash
# Linux (Debian/Ubuntu). El `sudo -v` es para que pida la clave antes del pipe.
sudo -v && curl https://rclone.org/install.sh | sudo bash

# macOS
brew install rclone

rclone version     # 1.60 o superior alcanza de sobra
```

No uses el `rclone` de `apt`: suele estar varias versiones atrás, y las de Drive
mejoraron bastante.

---

## 2. Un `client_id` propio — no te lo saltees

Esto es lo que más impacta y lo que casi todas las guías dejan como opcional.

rclone trae un `client_id` compartido por **todos** sus usuarios. Google aplica
la cuota por client_id, así que ese está permanentemente saturado. Con los
~190 GB de `.sra` y ~340 GB de BAMs de este proyecto, la diferencia es entre
horas y días, con errores `403 rateLimitExceeded` intermitentes en el medio.

Toma unos 5 minutos:

1. Andá a [console.cloud.google.com](https://console.cloud.google.com) con
   **la misma cuenta** (`seb.ugazm@gmail.com`).
2. **Crear un proyecto** — nombre cualquiera, p. ej. `tesis-rclone`.
3. **APIs y servicios → Biblioteca** → buscá *Google Drive API* → **Habilitar**.
4. **APIs y servicios → Pantalla de consentimiento de OAuth**:
   - Tipo: **Externo** (no hay Workspace).
   - Nombre de la app, tu email de soporte y tu email de contacto. Nada más.
   - En **Usuarios de prueba**, agregá `seb.ugazm@gmail.com`. Sin esto la
     autorización falla con *"app no verificada"*.
   - No hace falta publicar ni pedir verificación: con la cuenta como usuario
     de prueba alcanza. El token caduca cada 7 días en modo *Testing*; si te
     cansa, pasá la app a *En producción* (no requiere verificación mientras
     solo la uses vos con scopes de tu propia cuenta).
5. **Credenciales → Crear credenciales → ID de cliente de OAuth**:
   - Tipo de aplicación: **Aplicación de escritorio**. No *Aplicación web*:
     esa exige declarar un *redirect URI* y rclone usa un puerto local que
     cambia.
   - Copiá el **Client ID** y el **Client Secret** con el botón de copiar, no
     a mano. El Client ID **termina en `.apps.googleusercontent.com`**; si lo
     que copiaste no termina así, copiaste otra cosa (el nombre de la
     credencial, o el número de proyecto).

Los dos valores tardan **unos minutos** en propagarse del lado de Google. Si
autorizás en el mismo instante en que creaste la credencial, la primera vez
puede fallar y andar en el segundo intento.

---

## 3. `rclone config`

```bash
rclone config
```

Respuestas, en orden:

| pregunta | respuesta | por qué |
| :-- | :-- | :-- |
| `e/n/d/r/c/s/q>` | **`n`** | nuevo remoto |
| `name>` | **`gdrive-tesis`** | es el nombre que los scripts buscan |
| `Storage>` | **`drive`** | escribí la palabra, no el número: el número cambia entre versiones |
| `client_id>` | el del paso 2 | dejarlo vacío usa el compartido y saturado |
| `client_secret>` | el del paso 2 | |
| `scope>` | **`1`** (`drive`) | **ver abajo — acá es donde se rompe** |
| `service_account_file>` | *(vacío)* | |
| `Edit advanced config?` | **`n`** | |
| `Use web browser to automatically authenticate?` | **`y`** con navegador, **`n`** sin él | ver paso 4 |
| `Configure this as a Shared Drive (Team Drive)?` | **`n`** | es Mi unidad |
| `Keep this "gdrive-tesis" remote?` | **`y`** | |

### El `scope` es el error que no da error

`scope` ofrece varias opciones. Las dos que importan:

- **`1` — `drive`**: acceso completo. **Es la que va.**
- `3` — `drive.file`: rclone solo ve **los ficheros que él mismo creó**.

Las carpetas de `tesis/` (`80_sra/`, `10_bam/`, …) se crearon a mano en la web.
Con `drive.file` rclone **no las ve** — y no tira un error de permisos: el
remoto simplemente **lista vacío**. `drive_pull.sh sra prupe --go` baja 0
ficheros y sale con código 0, como si hubiera terminado bien.

`drive_check.sh` lo detecta y lo nombra.

### `root_folder_id`: no lo pongas

En *advanced config* se puede fijar `root_folder_id` al ID de `tesis/`
(`1E_Q6XLg4_RD01TFgfHbGUSQxx2ExtuVP`). **No lo hagas**, porque los scripts ya
anteponen `tesis/` por ruta: rclone terminaría buscando `tesis/tesis/80_sra/`,
que no existe — y otra vez, listando vacío en vez de fallar.

Si ya lo tenés puesto y no querés sacarlo, la salida es `export DRIVE_ROOT=''`.
`drive_check.sh` verifica que las dos cosas no coexistan.

---

## 4. Autorizar

### Con navegador en la misma máquina

`y` en *"Use web browser"*. rclone abre el navegador, elegís la cuenta, y en la
pantalla de *"Google no ha verificado esta aplicación"* → **Configuración
avanzada → Ir a (tu app)**. Es tu propia app, creada en el paso 2.

### Sin navegador (servidor, SSH)

Respondé **`n`**. rclone imprime un comando. En una máquina **con** navegador,
con rclone instalado:

```bash
rclone authorize "drive" "<CLIENT_ID>" "<CLIENT_SECRET>"
```

Autorizás ahí, copiás el bloque de token que imprime, y lo pegás en la máquina
sin navegador. No hace falta que la segunda máquina tenga la config: es un
one-shot.

---

## 4b. Cuando la autorización falla

El fallo de OAuth **muere en la pestaña del navegador y no vuelve a rclone**.
`rclone config` sigue como si nada y te ofrece `Keep this remote? y` — si decís
que sí, queda un remoto **sin token**, y a partir de ahí todos los comandos
fallan sin decir de dónde viene. Por eso `drive_check.sh` mira el token: es el
rastro que deja este problema.

### `Error 401: invalid_client` / *"The OAuth client was not found"*

Google no pudo resolver el `client_id` que rclone le mandó. **No es un problema
de permisos ni de la cuenta**: es que esa credencial, tal como llegó, no existe.
Cinco causas, en orden de frecuencia:

1. **El `client_id` está cortado o no es el `client_id`.** Tiene que terminar en
   `.apps.googleusercontent.com`. Un valor como `1234567890-a1b2c3` a secas da
   exactamente este error.
2. **Un espacio o un salto de línea pegados al copiar.** No se ve en la consola
   de rclone y Google lo manda igual.
3. **Pegaste el Client Secret en el campo `client_id`**, o al revés.
4. **La credencial está en otro proyecto de Google Cloud** que el de la Drive
   API, o se creó con otra cuenta de Google (fijate el selector de proyecto
   arriba a la izquierda en la consola).
5. **La credencial se borró**, o todavía no propagó (esperá un par de minutos y
   reintentá).

Para ver qué tiene guardado rclone:

```bash
rclone config dump | python3 -m json.tool | grep -A1 client_id
```

Compará carácter por carácter contra **Credenciales → tu cliente OAuth** en la
consola. Para corregirlo sin rehacer el remoto:

```bash
rclone config
# e -> gdrive-tesis -> pegá de nuevo client_id y client_secret
# -> al final dice "Already have a token - refresh?" -> y
```

### Otros dos que se confunden con este

| lo que dice el navegador | qué es |
| :-- | :-- |
| `Error 400: redirect_uri_mismatch` | creaste el cliente como **Aplicación web** en vez de **Aplicación de escritorio** |
| `Acceso bloqueado: ... no completó el proceso de verificación` | `seb.ugazm@gmail.com` no está en **Usuarios de prueba** de la pantalla de consentimiento (paso 2.4) |

Ninguno de los dos se arregla desde rclone: son de la consola de Google.

### Y el que aparece una semana después

Con la app en modo *Testing*, el refresh token **caduca a los 7 días** y el
próximo comando falla con `Token has been expired or revoked`. Se re-autoriza
con `rclone config reconnect gdrive-tesis:`. Para que no vuelva a pasar, pasá la
app a **En producción** en la pantalla de consentimiento — no requiere
verificación mientras la uses solo vos sobre tu propia cuenta.

---

## 5. Verificar

```bash
./scripts/drive_check.sh
```

Chequea nueve cosas: rclone instalado, el remoto existe y es `drive`, la
credencial OAuth (`client_id` bien formado y token presente), el `scope`, que
`root_folder_id` no choque con `DRIVE_ROOT`, la cuenta, que la raíz se vea, que
estén las 8 carpetas del mapa de fases, y el espacio libre.

Después, un dry-run de verdad — no baja nada:

```bash
./scripts/drive_pull.sh sra prupe
```

Tiene que listar ficheros. Si lista 0 y `drive_check.sh` dio verde, mirá que
`80_sra/prupe/` exista en Drive.

---

## 6. Dónde queda la credencial

```bash
rclone config file      # imprime la ruta
```

Normalmente `~/.config/rclone/rclone.conf`. Contiene un **refresh token**: quien
lo tenga entra a tu Drive completo hasta que lo revoques.

- Está en `.gitignore` (`rclone.conf` y `.rclone.conf`) — **no lo commitees**.
- `chmod 600 ~/.config/rclone/rclone.conf`.
- Para revocar: [myaccount.google.com/permissions](https://myaccount.google.com/permissions).
- Se puede cifrar con `rclone config` → `s` (*Set configuration password*).
  Pide la clave en cada corrida, así que no conviene para procesos largos.

---

## 7. Qué flags usan los scripts, y por qué

No hace falta que los pongas: están en el código. Pero conviene saber qué hacen.

| flag | dónde | por qué |
| :-- | :-- | :-- |
| `--checksum` | push, pull, purge | Compara hash, no tamaño+fecha. Un BAM re-generado con otro orden de `sort` puede pesar igual y ser distinto |
| `--transfers 4` | push, pull | Más paralelismo con Drive da más `403` que velocidad |
| `--drive-chunk-size 64M` | push | Menos llamadas por fichero grande. Cuesta 64 MB de RAM por transferencia: con `--transfers 4` son 256 MB |
| `--progress` | push, pull | |
| `--dry-run` | sin `--go` | **El default es dry-run.** Mover cientos de GB no puede ser lo que pasa si te equivocás de comando |
| `check --one-way` | `purge` | Confirma que Drive tiene lo local **antes** de borrar |

### Si te topás con límites

Google corta a ~750 GB/día de subida por cuenta. Los ~340 GB de BAMs entran, pero
si alguna vez chocás:

```bash
rclone ... --tpslimit 10 --retries 10 --low-level-retries 20
```

---

## 8. El primer uso real

```bash
./scripts/drive_check.sh                      # que la config esté bien
./scripts/drive_pull.sh sra prupe             # dry-run: qué bajaría
./scripts/drive_pull.sh sra prupe --go        # bajar de verdad
./scripts/trim.sh plan prupe                  # ¿los encuentra el recorte?
./scripts/trim.sh correr prupe                # recortar
# ... alinear ...
./scripts/drive_push.sh bam prupe --go        # subir los BAMs
./scripts/drive_pull.sh purge sra prupe --go  # liberar disco
```

`drive_pull.sh` deja los `.sra` en el mismo directorio donde `trim.sh` los
busca — eso lo garantiza `ruta_local()` en `_drive_lib.sh` y lo verifica
`tests/test_rutas.sh`. No siempre fue así: hubo un momento en que tres scripts
apuntaban a tres lugares distintos.

**`purge` se niega a borrar** si Drive no tiene todo lo local, salvo para `sra`
y `genomas`, donde Drive es la fuente. Para `bam` la copia local es lo que
produjo el alineamiento: si todavía no se subió, borrarla la pierde.
